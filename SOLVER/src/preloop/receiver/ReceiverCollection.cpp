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

ReceiverCollection::ReceiverCollection(const std::string &fileRec, bool geographic, 
    double srcLat, double srcLon, double srcDep):
mInputFile(fileRec), mGeographic(geographic) {
    std::vector<std::string> name, network;
    std::vector<double> theta, phi, depth;
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
                if (strs.size() < 5 || strs.size() > 6) {
                    continue;
                }
                name.push_back(strs[0]);
                network.push_back(strs[1]);
                theta.push_back(boost::lexical_cast<double>(strs[2]));
                phi.push_back(boost::lexical_cast<double>(strs[3]));
                depth.push_back(boost::lexical_cast<double>(strs[strs.size() - 1]));
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
    
    // create receivers
    mWidthName = -1;
    mWidthNetwork = -1;
    std::vector<std::string> recKeys;
    for (int i = 0; i < name.size(); i++) {
        mReceivers.push_back(new Receiver(name[i], network[i], 
            theta[i], phi[i], geographic, depth[i], srcLat, srcLon, srcDep));
        mWidthName = std::max(mWidthName, (int)name[i].length());
        mWidthNetwork = std::max(mWidthNetwork, (int)network[i].length());
        // check duplicated
        std::string key = network[i] + "_" + name[i];
        if (std::find(recKeys.begin(), recKeys.end(), key) != recKeys.end()) {
            throw std::runtime_error("ReceiverCollection::ReceiverCollection || "
                "Duplicated station keys (network_name) found in station data file " + mInputFile + " || "
                "Name = " + name[i] + "; Network = " + network[i]);
        }
        recKeys.push_back(key);
    }
}

ReceiverCollection::~ReceiverCollection() {
    for (const auto &rec: mReceivers) {
        delete rec;
    }
}

void ReceiverCollection::release(Domain &domain, const Mesh &mesh) {
    // locate receivers
    MultilevelTimer::begin("Locate Receivers", 2);
    std::vector<int> recRank(mReceivers.size(), XMPI::nproc());
    std::vector<int> recETag(mReceivers.size(), -1);
    std::vector<RDMatPP> recInterpFact(mReceivers.size(), RDMatPP::Zero());
    for (int irec = 0; irec < mReceivers.size(); irec++) {
        bool found = mReceivers[irec]->locate(mesh, recETag[irec], recInterpFact[irec]);
        if (found) {
            recRank[irec] = XMPI::rank();
        }
    }
    MultilevelTimer::end("Locate Receivers", 2);
    
    // release to domain
    MultilevelTimer::begin("Release to Domain", 2);
    PointwiseRecorder *recorderPW = new PointwiseRecorder(
        mTotalRecordSteps, mRecordInterval, mBufferSize, mENZ);
    for (int irec = 0; irec < mReceivers.size(); irec++) {
        int recRankMin = XMPI::min(recRank[irec]);
        if (recRankMin == XMPI::nproc()) {
            throw std::runtime_error("ReceiverCollection::release || Error locating receiver || " 
                "Name = " + mReceivers[irec]->getName() + "; "
                "Network = " + mReceivers[irec]->getNetwork());
        }
        if (recRankMin == XMPI::rank()) {
            mReceivers[irec]->release(*recorderPW, 
                domain, recETag[irec], recInterpFact[irec]);
        }
    }
    
    // IO
    if (mAscii) {
        recorderPW->addIO(new PointwiseIOAscii());
    }
    if (mNetCDF) {
        recorderPW->addIO(new PointwiseIONetCDF());
    }
    if (mASDF) {
        throw std::runtime_error("ReceiverCollection::release || ASDF not implemented.");
    }
    
    // add recorder to domain
    domain.setPointwiseRecorder(recorderPW);
    MultilevelTimer::begin("Release to Domain", 2);
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
    ss << "========================= Receivers ========================\n" << std::endl;
    return ss.str();
}

void ReceiverCollection::buildInparam(ReceiverCollection *&rec, const Parameters &par, 
    double srcLat, double srcLon, double srcDep, int totalStepsSTF, int verbose) {
    if (rec) {
        delete rec;
    }
    
    // create from file
    std::string recFile = Parameters::sInputDirectory + "/" 
        + par.getValue<std::string>("OUT_STATIONS_FILE");
    std::string recSys = par.getValue<std::string>("OUT_STATIONS_SYSTEM");
    bool geographic;
    if (boost::iequals(recSys, "source_centered")) {
        geographic = false;
    } else if (boost::iequals(recSys, "geographic")) {
        geographic = true;
    } else {
        throw std::runtime_error("ReceiverCollection::buildInparam || "
            "Invalid parameter, keyword = OUT_STATIONS_SYSTEM.");
    }
    rec = new ReceiverCollection(recFile, geographic, srcLat, srcLon, srcDep); 
    
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
        rec->mENZ = false;
    } else if (boost::iequals(strcomp, "ENZ")) {
        rec->mENZ = true;
    } else {
        throw std::runtime_error("ReceiverCollection::buildInparam || "
            "Invalid parameter, keyword = OUT_STATIONS_COMPONENTS.");
    }
    
    // IO
    int numFmt = par.getSize("OUT_STATIONS_FORMAT");
    for (int i = 0; i < numFmt; i++) {
        std::string strfmt = par.getValue<std::string>("OUT_STATIONS_FORMAT", i); 
        if (boost::iequals(strfmt, "ascii")) {
            rec->mAscii = true;
        } else if (boost::iequals(strfmt, "netcdf")) {
            rec->mNetCDF = true;
        } else if (boost::iequals(strfmt, "asdf")) {
            rec->mASDF = true;
        } else {
            throw std::runtime_error("ReceiverCollection::buildInparam || "
                "Invalid parameter, keyword = OUT_STATIONS_FORMAT.");
        }
    }
    
    if (verbose) {
        XMPI::cout << rec->verbose();
    }
}

