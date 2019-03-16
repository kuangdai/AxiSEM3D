// SurfaceRecorder.cpp
// created by Kuangdai on 27-Nov-2017 
// recorder for surface wavefield

#include "SurfaceRecorder.h"
#include "SurfaceIO.h"
#include "SurfaceInfo.h"
#include "XMPI.h"

SurfaceRecorder::SurfaceRecorder(int totalRecordSteps, int recordInterval, 
    int bufferSize, double srcLat, double srcLon, double srcDep, bool assemble): 
mTotalRecordSteps(totalRecordSteps),
mRecordInterval(recordInterval), mBufferSize(bufferSize),
mSrcLat(srcLat), mSrcLon(srcLon), mSrcDep(srcDep) {
    mBufferLine = 0;
    mIO = new SurfaceIO(assemble);
}

SurfaceRecorder::~SurfaceRecorder() {
    delete mIO;
}

void SurfaceRecorder::addElement(Element *ele, int surfSide) {
    mSurfaceInfo.push_back(SurfaceInfo(ele, surfSide));
}

void SurfaceRecorder::initialize() {
    // global tag
    int numSurfEle = mSurfaceInfo.size();
    std::vector<int> numSurfEleAll;
    XMPI::gather(numSurfEle, numSurfEleAll, true);
    int startGlobalTag = 0;
    for (int iproc = 0; iproc < XMPI::rank(); iproc++) {
        startGlobalTag += numSurfEleAll[iproc];
    }
    // theta and tag
    std::vector<double> minTheta;
    std::vector<int> unsortedTag;
    for (int iele = 0; iele < numSurfEle; iele++) {
        minTheta.push_back(mSurfaceInfo[iele].getThetaMin());
        unsortedTag.push_back(startGlobalTag + iele);
    }
    // gather
    std::vector<std::vector<double>> minThetaAll;
    std::vector<std::vector<int>> unsortedTagAll;
    XMPI::gather(minTheta, minThetaAll, MPI_DOUBLE, true);
    XMPI::gather(unsortedTag, unsortedTagAll, MPI_INT, true);
    // flatten
    std::vector<std::pair<double, int>> minThetaTagAll;
    for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
        for (int iele = 0; iele < numSurfEleAll[iproc]; iele++) {
            minThetaTagAll.push_back(std::make_pair(minThetaAll[iproc][iele], 
                unsortedTagAll[iproc][iele]));
        }
    }
    // sort
    std::sort(minThetaTagAll.begin(), minThetaTagAll.end());
    // assigned sorted tag
    for (int gtag = 0; gtag < minThetaTagAll.size(); gtag++) {
        int ltag = minThetaTagAll[gtag].second - startGlobalTag;
        if (ltag >= 0 && ltag < numSurfEle) {
            mSurfaceInfo[ltag].setGlobalTag(gtag);
        }
    }
    
    // buffer
    mBufferTime = RDColX::Zero(mBufferSize);
    for (int iele = 0; iele < numSurfEle; iele++) {
        CMatXX_RM buf;
        mSurfaceInfo[iele].initBuffer(mBufferSize, buf);
        mBufferDisp.push_back(buf);
    }    
    
    // IO
    mIO->initialize(mTotalRecordSteps, mBufferSize, mSurfaceInfo, 
        mSrcLat, mSrcLon, mSrcDep);
}

void SurfaceRecorder::finalize() {
    mIO->finalize();
}

void SurfaceRecorder::record(int tstep, double t) {
    if (tstep % mRecordInterval != 0) {
        return;
    }
    
    // time
    mBufferTime(mBufferLine) = t;
    
    // get disp
    for (int iele = 0; iele < mSurfaceInfo.size(); iele++) {
        // compute from element
        mSurfaceInfo[iele].feedBuffer(mBufferLine, mBufferDisp[iele]);
    }
    
    // increment buffer line
    mBufferLine++;
    
    // dump and clear buffer
    if (mBufferLine == mBufferSize) {
        dumpToFile();
    }
}

void SurfaceRecorder::dumpToFile() {
    mIO->dumpToFile(mBufferDisp, mBufferTime, mBufferLine);
    mBufferLine = 0;
}

