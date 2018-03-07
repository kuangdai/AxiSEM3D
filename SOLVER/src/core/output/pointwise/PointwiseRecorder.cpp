// PointwiseRecorder.h
// created by Kuangdai on 1-Jun-2017 
// recorder for point-wise receivers

#include "PointwiseRecorder.h"
#include "Element.h"
#include "PointwiseIO.h"

PointwiseRecorder::PointwiseRecorder(int totalRecordSteps, int recordInterval, 
    int bufferSize, const std::string &components, 
    double srcLat, double srcLon, double srcDep): 
mTotalRecordSteps(totalRecordSteps),
mRecordInterval(recordInterval), mBufferSize(bufferSize), mComponents(components),
mSrcLat(srcLat), mSrcLon(srcLon), mSrcDep(srcDep) {
    mBufferLine = 0;
}

PointwiseRecorder::~PointwiseRecorder() {
    for (const auto &io: mIOs) {
        delete io;
    }
}

void PointwiseRecorder::addReceiver(const std::string &name, const std::string &network, 
    double phi, const RDMatPP &weights, const Element *ele, double theta, double baz,
    double lat, double lon, double dep, bool dumpStrain) {
    mPointwiseInfo.push_back(PointwiseInfo(name, network, phi, weights, ele, theta, baz, 
        lat, lon, dep, dumpStrain));
}

void PointwiseRecorder::initialize() {
    int numRec = mPointwiseInfo.size();
    mBufferDisp = RMatXX_RM::Zero(mBufferSize, numRec * 3);
    mBufferTime = RColX::Zero(mBufferSize);
    int numStrainRec = 0;
    for (const auto &rec: mPointwiseInfo) {
        if (rec.mDumpStrain) {
            numStrainRec++;
        }
    }
    mBufferStrain = RMatXX_RM::Zero(mBufferSize, numStrainRec * 6);
    for (const auto &io: mIOs) {
        io->initialize(mTotalRecordSteps, mBufferSize, mComponents, mPointwiseInfo,
            mSrcLat, mSrcLon, mSrcDep);
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
        if (mComponents != "SPZ") {
            // transform
            Real cost = cos(mPointwiseInfo[irec].mTheta);
            Real sint = sin(mPointwiseInfo[irec].mTheta);
            Real ur = gm(0) * sint + gm(2) * cost; 
            Real ut = gm(0) * cost - gm(2) * sint;
            if (mComponents == "ENZ") {
                Real cosbaz = cos(mPointwiseInfo[irec].mBAz);
                Real sinbaz = sin(mPointwiseInfo[irec].mBAz);
                gm(0) = -ut * sinbaz + gm(1) * cosbaz;
                gm(1) = -ut * cosbaz - gm(1) * sinbaz;
                gm(2) = ur;
            } else { 
                // RTZ
                gm(0) = ut;
                gm(2) = ur;
            }
        }    
        // write to buffer
        mBufferDisp.block(mBufferLine, irec * 3, 1, 3) = gm;
    }
    
    // get disp
    static RRow6 strain;
    int istrain = 0;
    for (int irec = 0; irec < mPointwiseInfo.size(); irec++) {
        if (mPointwiseInfo[irec].mDumpStrain) {
            // compute from element
            mPointwiseInfo[irec].mElement->computeStrain(mPointwiseInfo[irec].mPhi, 
                mPointwiseInfo[irec].mWeights, strain);
            // write to buffer
            mBufferStrain.block(mBufferLine, istrain * 6, 1, 6) = strain;
            istrain++;
        }
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
        io->dumpToFile(mBufferDisp, mBufferStrain, mBufferTime, mBufferLine);
    }
    mBufferLine = 0;
}

