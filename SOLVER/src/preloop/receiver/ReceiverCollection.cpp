// Receiver.cpp
// created by Kuangdai on 1-Jun-2016 
// receiver collections

#include "ReceiverCollection.h"
#include "Receiver.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include "XMPI.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "Parameters.h"
#include "Domain.h"
#include "MultilevelTimer.h"
#include "PointwiseRecorder.h"
#include "PointwiseIOAscii.h"
#include "PointwiseIONetCDF.h"
#include "SurfaceRecorder.h"
#include "Mesh.h"
#include "Quad.h"

ReceiverCollection::ReceiverCollection(const std::string &fileRec, bool geographic, 
    double srcLat, double srcLon, double srcDep, int duplicated, bool saveSurf):
mInputFile(fileRec), mGeographic(geographic), mSaveWholeSurface(saveSurf),
mSrcLat(srcLat), mSrcLon(srcLon), mSrcDep(srcDep) {
    std::vector<std::string> name, network;
    std::vector<double> theta, phi, depth;
    std::vector<int> dumpStrain;
    std::vector<int> dumpCurl;
    if (!boost::iequals(fileRec, "none")) {
        mInputFile = Parameters::sInputDirectory + "/" + mInputFile;
        if (XMPI::root()) {
            std::fstream fs(mInputFile, std::fstream::in);
            if (!fs) {
                throw std::runtime_error("ReceiverCollection::ReceiverCollection || "
                    "Error opening station data file " + mInputFile + ".");
            }
            std::string line;
            while (getline(fs, line)) {
                try {
                    std::vector<std::string> strs = Parameters::splitString(line, "\t ");
                    int npar = strs.size();
                    if (npar < 6) {
                        continue;
                    }
                    double theta_try = boost::lexical_cast<double>(strs[2]);
                    double phi_try = boost::lexical_cast<double>(strs[3]);
                    // the 4th column (elevation) is ignored
                    double depth_try = boost::lexical_cast<double>(strs[5]);
                    int strain_try = 0;
                    int curl_try = 0;
                    if (npar > 6) {
                        for (int ipar = 6; ipar < npar; ipar++) {
                            if (boost::iequals(strs[ipar], "dump_strain")) {
                                strain_try = 1;
                            }
                            if (boost::iequals(strs[ipar], "dump_curl")) {
                                curl_try = 1;
                            }
                        }
                    }
                    name.push_back(strs[0]);
                    network.push_back(strs[1]);
                    theta.push_back(theta_try);
                    phi.push_back(phi_try);
                    depth.push_back(depth_try);
                    dumpStrain.push_back(strain_try);
                    dumpCurl.push_back(curl_try);
                } catch(std::exception) {
                    // simply ignore invalid lines
                    continue;
                }
            }
            fs.close();
        }
        XMPI::bcast(name);
        XMPI::bcast(network);
        XMPI::bcast(theta);
        XMPI::bcast(phi);
        XMPI::bcast(depth);
        XMPI::bcast(dumpStrain);
        XMPI::bcast(dumpCurl);
    }
    
    // create receivers
    mWidthName = -1;
    mWidthNetwork = -1;
    for (int i = 0; i < name.size(); i++) {
        // check duplicated
        if (duplicated != 0) {
            std::vector<std::string> recKeys;
            std::string key = network[i] + "." + name[i];
            if (std::find(recKeys.begin(), recKeys.end(), key) != recKeys.end()) {
                if (duplicated == 1) {
                    // rename
                    int append = 0;
                    std::string nameOriginal = name[i];
                    while (std::find(recKeys.begin(), recKeys.end(), key) != recKeys.end()) {
                        name[i] = nameOriginal + "__DUPLICATED" + boost::lexical_cast<std::string>(++append);
                        key = network[i] + "." + name[i];
                    }
                } else {
                    // error
                    throw std::runtime_error("ReceiverCollection::ReceiverCollection || "
                        "Duplicated station keys (network_name) found in station data file " + mInputFile + " || "
                        "Name = " + name[i] + "; Network = " + network[i]);
                }
            }
            recKeys.push_back(key);
        }
        // add receiver
        mReceivers.push_back(new Receiver(name[i], network[i], 
            theta[i], phi[i], geographic, depth[i], (bool)dumpStrain[i], (bool)dumpCurl[i], 
            srcLat, srcLon, srcDep));
        mWidthName = std::max(mWidthName, (int)name[i].length());
        mWidthNetwork = std::max(mWidthNetwork, (int)network[i].length());
    }
}

ReceiverCollection::~ReceiverCollection() {
    for (const auto &rec: mReceivers) {
        delete rec;
    }
}

void ReceiverCollection::release(Domain &domain, const Mesh &mesh, bool depthInRef) {
    // locate receivers
    MultilevelTimer::begin("Locate Receivers", 2);
    std::vector<int> recRank(mReceivers.size(), XMPI::nproc());
    std::vector<int> recETag(mReceivers.size(), -1);
    std::vector<int> recQTag(mReceivers.size(), -1);
    // std::vector<RDMatPP> recInterpFact(mReceivers.size(), RDMatPP::Zero());
    for (int irec = 0; irec < mReceivers.size(); irec++) {
        bool found = mReceivers[irec]->locate(mesh, recETag[irec], recQTag[irec], depthInRef);
        if (found) {
            recRank[irec] = XMPI::rank();
        }
    }
    MultilevelTimer::end("Locate Receivers", 2);
    
    // release to domain
    MultilevelTimer::begin("Release to Domain", 2);
    PointwiseRecorder *recorderPW = new PointwiseRecorder(
        mTotalRecordSteps, mRecordInterval, mBufferSize, mComponents,
        mSrcLat, mSrcLon, mSrcDep);
    
    MultilevelTimer::begin("Find Min Rank", 3);    
    std::vector<int> recRankMinG(recRank);
    XMPI::min(recRank, recRankMinG);    
    MultilevelTimer::end("Find Min Rank", 3);
    
    MultilevelTimer::begin("Release receivers", 3);
    for (int irec = 0; irec < mReceivers.size(); irec++) {
        int recRankMin = recRankMinG[irec];
        if (recRankMin == XMPI::nproc()) {
            throw std::runtime_error("ReceiverCollection::release || Error locating receiver || " 
                "Name = " + mReceivers[irec]->getName() + "; "
                "Network = " + mReceivers[irec]->getNetwork());
        }
        if (recRankMin == XMPI::rank()) {
            RDMatPP interpFact;
            mReceivers[irec]->computeInterpFact(mesh, recQTag[irec], interpFact, depthInRef);
            mReceivers[irec]->release(*recorderPW, 
                domain, recETag[irec], interpFact);
        }
    }
    MultilevelTimer::end("Release receivers", 3);
    
    // IO
    for (const auto &io: mPointwiseIO) {
        recorderPW->addIO(io);
    }
    
    // add recorder to domain
    domain.setPointwiseRecorder(recorderPW);
    
    // whole surface
    if (mSaveWholeSurface) {
        MultilevelTimer::begin("Whole Surface", 3);
        SurfaceRecorder *recorderSF = new SurfaceRecorder(mTotalRecordSteps, 
            mRecordInterval, mBufferSize, 
            mSrcLat, mSrcLon, mSrcDep, mAssemble);
        for (int iloc = 0; iloc < mesh.getNumQuads(); iloc++) {
            const Quad *quad = mesh.getQuad(iloc);
            if (quad->onSurface()) {
                Element *ele = domain.getElement(quad->getElementTag());
                recorderSF->addElement(ele, quad->getSurfSide());
            }
        }
        domain.setSurfaceRecorder(recorderSF);
        MultilevelTimer::end("Whole Surface", 3);
    }
    MultilevelTimer::end("Release to Domain", 2);
}

std::string ReceiverCollection::verbose() const {
    std::stringstream ss;
    ss << "\n========================= Receivers ========================" << std::endl;
    ss << "  Number of Receivers   =   " << mReceivers.size() << std::endl;
    ss << "  Coordinate System     =   " << (mGeographic ? "Geographic" : "Source-centered") << std::endl;
    if (mReceivers.size() > 0) {
        ss << "  Receiver List: " << std::endl;
        ss << "    " << mReceivers[0]->verbose(mGeographic, mWidthName, mWidthNetwork) << std::endl;
    }
    if (mReceivers.size() > 1) {
        ss << "    " << std::setw(mWidthName) << "..." << std::endl;
        ss << "    " << mReceivers[mReceivers.size() - 1]->verbose(mGeographic, mWidthName, mWidthNetwork) << std::endl;
    }
    if (mSaveWholeSurface) {
        ss << "  * Wavefield on the whole surface will be saved." << std::endl;
    }
    ss << "========================= Receivers ========================\n" << std::endl;
    return ss.str();
}

void ReceiverCollection::buildInparam(ReceiverCollection *&rec, const Parameters &par, 
    double srcLat, double srcLon, double srcDep, int totalStepsSTF, int verbose) {
    if (rec) {
        delete rec;
    }
    
    // create from file
    std::string recFile = par.getValue<std::string>("OUT_STATIONS_FILE");
    std::string recSys = par.getValue<std::string>("OUT_STATIONS_SYSTEM");
    bool geographic;
    if (boost::iequals(recSys, "source-centered")) {
        geographic = false;
    } else if (boost::iequals(recSys, "geographic")) {
        geographic = true;
    } else {
        throw std::runtime_error("ReceiverCollection::buildInparam || "
            "Invalid parameter, keyword = OUT_STATIONS_SYSTEM.");
    }
    std::string dupOption = par.getValue<std::string>("OUT_STATIONS_DUPLICATED");
    int duplicated = 0;
    if (boost::iequals(dupOption, "ignore")) {
        duplicated = 0;
    } else if (boost::iequals(dupOption, "rename")) {
        duplicated = 1;
    } else if (boost::iequals(dupOption, "error")) {
        duplicated = 2;
    } else {
        throw std::runtime_error("ReceiverCollection::buildInparam || "
            "Invalid parameter, keyword = OUT_STATIONS_DUPLICATED.");
    }
    bool saveSurf = par.getValue<bool>("OUT_STATIONS_WHOLE_SURFACE");
    rec = new ReceiverCollection(recFile, geographic, srcLat, srcLon, srcDep, 
        duplicated, saveSurf); 
    
    // options 
    rec->mRecordInterval = par.getValue<int>("OUT_STATIONS_RECORD_INTERVAL");
    if (rec->mRecordInterval <= 0) {
        rec->mRecordInterval = 1;
    }
    rec->mTotalRecordSteps = totalStepsSTF / rec->mRecordInterval;
    if (totalStepsSTF % rec->mRecordInterval > 0) {
        rec->mTotalRecordSteps += 1;
    }
    rec->mBufferSize = par.getValue<int>("OUT_STATIONS_DUMP_INTERVAL");
    if (rec->mBufferSize <= 0) {
        rec->mBufferSize = 1000;
    }
    if (rec->mBufferSize > rec->mTotalRecordSteps) {
        rec->mBufferSize = rec->mTotalRecordSteps;
    }
    std::string strcomp = par.getValue<std::string>("OUT_STATIONS_COMPONENTS");
    if (boost::iequals(strcomp, "RTZ")) {
        rec->mComponents = "RTZ";
    } else if (boost::iequals(strcomp, "ENZ")) {
        rec->mComponents = "ENZ";
    } else if (boost::iequals(strcomp, "SPZ")) {
        rec->mComponents = "SPZ";
    } else {
        throw std::runtime_error("ReceiverCollection::buildInparam || "
            "Invalid parameter, keyword = OUT_STATIONS_COMPONENTS.");
    }
    
    // IO
    int numFmt = par.getSize("OUT_STATIONS_FORMAT");
    // use bool first to avoid duplicated IO
    bool ascii = false, netcdf = false, netcdf_no_assemble = false; 
    for (int i = 0; i < numFmt; i++) {
        std::string strfmt = par.getValue<std::string>("OUT_STATIONS_FORMAT", i); 
        if (boost::iequals(strfmt, "ascii")) {
            ascii = true;
        } else if (boost::iequals(strfmt, "netcdf")) {
            netcdf = true;
        } else if (boost::iequals(strfmt, "netcdf_no_assemble")) {
            netcdf_no_assemble = true;
        } else {
            throw std::runtime_error("ReceiverCollection::buildInparam || "
                "Invalid parameter, keyword = OUT_STATIONS_FORMAT.");
        }
    }
    // create IO
    if (ascii) {
        rec->mPointwiseIO.push_back(new PointwiseIOAscii());
    }
    
    // check netcdf options
    if (netcdf && netcdf_no_assemble) {
        throw std::runtime_error("ReceiverCollection::buildInparam || "
            "Invalid parameter OUT_STATIONS_FORMAT, ||"
            "netcdf and netcdf_no_assemble cannot coexist.");
    }
    #ifdef _USE_PARALLEL_NETCDF
        if (netcdf_no_assemble) {
            throw std::runtime_error("ReceiverCollection::buildInparam || "
                "Invalid parameter OUT_STATIONS_FORMAT, ||"
                "netcdf_no_assemble is not compatible with parallel NetCDF.");
        }
    #endif
    
    if (netcdf) {
        rec->mPointwiseIO.push_back(new PointwiseIONetCDF(true));
    }
    if (netcdf_no_assemble) {
        rec->mPointwiseIO.push_back(new PointwiseIONetCDF(false));
        rec->mAssemble = false;
    }
    
    if (verbose) {
        XMPI::cout << rec->verbose();
    }
}

