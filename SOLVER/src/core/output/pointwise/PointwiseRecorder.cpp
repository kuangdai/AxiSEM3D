// PointwiseRecorder.h
// created by Kuangdai on 1-Jun-2017 
// recorder for point-wise receivers

#include "PointwiseRecorder.h"
#include "Element.h"
#include "PointwiseIO.h"

PointwiseRecorder::PointwiseRecorder(int totalRecordSteps, int recordInterval, 
    int bufferSize, bool ENZ): mTotalRecordSteps(totalRecordSteps),
mRecordInterval(recordInterval), mBufferSize(bufferSize), mENZ(ENZ) {
    mBufferLine = 0;
}

PointwiseRecorder::~PointwiseRecorder() {
    for (const auto &io: mIOs) {
        delete io;
    }
}

void PointwiseRecorder::addReceiver(const std::string &name, const std::string &network, 
    double phi, const RDMatPP &weights, const Element *ele, double theta, double baz,
    double lat, double lon, double dep) {
    mPointwiseInfo.push_back(PointwiseInfo(name, network, phi, weights, ele, theta, baz, lat, lon, dep));
}

void PointwiseRecorder::initialize() {
    int numRec = mPointwiseInfo.size();
    mBufferDisp = RMatXX_RM::Zero(mBufferSize, numRec * 3);
    mBufferTime = RColX::Zero(mBufferSize);
    std::vector<std::string> names;
    std::vector<std::string> networks;
    for (const auto &rec: mPointwiseInfo) {
        names.push_back(rec.mName);
        networks.push_back(rec.mNetwork);
    }
    for (const auto &io: mIOs) {
        io->initialize(mTotalRecordSteps, mBufferSize, mENZ, mPointwiseInfo);
    }
}

void PointwiseRecorder::finalize() {
    for (const auto &io: mIOs) {
        io->finalize();
    }
}

void PointwiseRecorder::record(int tstep, Real t) {
    if (tstep % mRecordInterval != 0) {
        return;
    }
    
    // time
    mBufferTime(mBufferLine) = t;
    
    // get disp
    static RRow3 gm;
    for (int irec = 0; irec < mPointwiseInfo.size(); irec++) {
        // compute from element
        mPointwiseInfo[irec].mElement->computeGroundMotion(mPointwiseInfo[irec].mPhi, 
            mPointwiseInfo[irec].mWeights, gm);
        // transform
        Real cost = cos(mPointwiseInfo[irec].mTheta);
        Real sint = sin(mPointwiseInfo[irec].mTheta);
        Real ur = gm(0) * sint + gm(2) * cost; 
        Real ut = gm(0) * cost - gm(2) * sint;
        if (mENZ) {
            Real cosbaz = cos(mPointwiseInfo[irec].mBAz);
            Real sinbaz = sin(mPointwiseInfo[irec].mBAz);
            gm(0) = -ut * sinbaz + gm(1) * cosbaz;
            gm(1) = -ut * cosbaz - gm(1) * sinbaz;
            gm(2) = ur;
        } else {
            gm(0) = ut;
            gm(2) = ur;
        }
        // write to buffer
        mBufferDisp.block(mBufferLine, irec * 3, 1, 3) = gm;
    }
    
    // increment buffer line
    mBufferLine++;
    
    // dump and clear buffer
    if (mBufferLine == mBufferSize) {
        dumpToFile();
    }
}

void PointwiseRecorder::dumpToFile() {
    for (const auto &io: mIOs) {
        io->dumpToFile(mBufferDisp, mBufferTime, mBufferLine);
    }
    mBufferLine = 0;
}

