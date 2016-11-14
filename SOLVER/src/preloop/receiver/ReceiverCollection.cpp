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

ReceiverCollection::ReceiverCollection(const std::string &fileRec, bool geographic, 
    double srcLat, double srcLon, double srcDep):
mInputFile(fileRec), mGeographic(geographic) {
    std::vector<std::string> name, network;
    std::vector<double> theta, phi, depth;
    if (XMPI::root()) {
        std::fstream fs(mInputFile, std::fstream::in);
        if (!fs) throw std::runtime_error("ReceiverCollection::ReceiverCollection || "
            "Error opening receiver data file " + mInputFile + ".");
        std::string line;
        while (getline(fs, line)) {
            try {
                std::vector<std::string> strs;
                boost::trim_if(line, boost::is_any_of("\t "));
                boost::split(strs, line, boost::is_any_of("\t "), boost::token_compress_on);
                if (strs.size() < 5 || strs.size() > 6) continue;
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
    for (int i = 0; i < name.size(); i++) {
        mReceivers.push_back(new Receiver(name[i], network[i], 
            theta[i], phi[i], geographic, depth[i], srcLat, srcLon, srcDep));
        mWidthName = std::max(mWidthName, (int)name[i].length());
        mWidthNetwork = std::max(mWidthNetwork, (int)network[i].length());
    }
}

ReceiverCollection::~ReceiverCollection() {
    for (const auto &rec: mReceivers)  delete rec;
}

void ReceiverCollection::release(Domain &domain, const Mesh &mesh) {
    std::vector<int> recRank(mReceivers.size(), XMPI::nproc());
    std::vector<int> recETag(mReceivers.size(), -1);
    std::vector<RDMatPP> recInterpFact(mReceivers.size(), RDMatPP::Zero());
    for (int irec = 0; irec < mReceivers.size(); irec++) {
        bool found = mReceivers[irec]->locate(mesh, recETag[irec], recInterpFact[irec]);
        if (found) recRank[irec] = XMPI::rank();
    }
    for (int irec = 0; irec < mReceivers.size(); irec++) {
        int recRankMin = XMPI::min(recRank[irec]);
        if (recRankMin == XMPI::nproc()) {
            throw std::runtime_error("ReceiverCollection::release || Error locating receiver " + 
                boost::lexical_cast<std::string>(irec));
        }
        if (recRankMin == XMPI::rank()) {
            mReceivers[irec]->release(domain, mesh, mRecordInterval, mComponent, 
                mOutputDir + "/stations", mBinary, mAppend, mBufferSize, 
                recETag[irec], recInterpFact[irec]);
        }
    }
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

void ReceiverCollection::buildInparam(ReceiverCollection *&rec, 
    const Parameters &par, double srcLat, double srcLon, double srcDep, int verbose) {
    if (rec) delete rec;
    
    // create from file
    std::string recFile = Parameters::sInputDirectory + "/STATIONS";
    std::string recSys = par.getValue<std::string>("OUT_STATIONS_SYSTEM");
    bool geographic;
    if (boost::iequals(recSys, "source_centered")) {
        geographic = false;
    } else if (boost::iequals(recSys, "geographic")) {
        geographic = true;
    } else {
        throw std::runtime_error("ReceiverCollection::buildInparam || Invalid parameter, keyword = OUT_STATIONS_SYSTEM.");
    }
    rec = new ReceiverCollection(recFile, geographic, srcLat, srcLon, srcDep); 
    
    // options 
    rec->mRecordInterval = par.getValue<int>("OUT_RECORD_INTERVAL");
    if (rec->mRecordInterval <= 0) rec->mRecordInterval = 1;
    std::string strcomp = par.getValue<std::string>("OUT_SEISMOGRAM_COMPONENTS");
    if (boost::iequals(strcomp, "RTZ")) {
        rec->mComponent = 0;
    } else if (boost::iequals(strcomp, "ENZ")) {
        rec->mComponent = 1;
    } else if (boost::iequals(strcomp, "SPZ")) {
        rec->mComponent = 2;
    } else {
        throw std::runtime_error("ReceiverCollection::buildInparam || Invalid parameter, keyword = OUT_SEISMOGRAM_COMPONENTS.");
    }
    rec->mOutputDir = Parameters::sOutputDirectory; 
    std::string strfmt = par.getValue<std::string>("OUT_SEISMOGRAM_FORMAT"); 
    if (boost::iequals(strfmt, "ascii")) {
        rec->mBinary = false;
    } else if (boost::iequals(strfmt, "binary")) {
        rec->mBinary = true;
    } else {
        throw std::runtime_error("ReceiverCollection::buildInparam || Invalid parameter, keyword = OUT_SEISMOGRAM_FORMAT.");
    }
    rec->mAppend = false;
    rec->mBufferSize = par.getValue<int>("OUT_DUMP_INTERVAL");
    if (rec->mBufferSize <= 0) rec->mBufferSize = 100;
    
    if (verbose) XMPI::cout << rec->verbose();
}

