// Mesh.cpp
// created by Kuangdai on 5-May-2016 
// AxiSEM3D Mesh

#include "Mesh.h"
#include "ExodusModel.h"
#include "Parameters.h"

#include "Connectivity.h"
#include "DualGraph.h"
#include "GLLPoint.h"
#include "Quad.h"
#include "XMPI.h"
#include "Domain.h"
#include "SolidElement.h"
#include "FluidElement.h"
#include "Point.h"
#include "XMath.h"
#include "Geometric3D.h"
#include "NuWisdom.h"

#include "XTimer.h"
#include "SlicePlot.h"
#include <fstream>

Mesh::~Mesh() {
    destroy(); // local build
    delete mDDPar;
    delete mLearnPar;
    for (const auto &sp: mSlicePlots) delete sp;    
}

Mesh::Mesh(const ExodusModel *exModel, const NrField *nrf, 
    double srcLat, double srcLon, double srcDep, const Parameters &par):
mExModel(exModel), mNrField(nrf), mSrcLat(srcLat), mSrcLon(srcLon), mSrcDep(srcDep) {
    mAttBuilder = 0;
    mMsgInfo = 0;
    mOceanLoad3D = 0;
    mDDPar = new DDParameters(par);
    mLearnPar = new LearnParameters(par);
    
    // 2D mode
    mUse2D = par.getValue<bool>("MODEL_2D_MODE");
    mPhi2D = -1.;
    if (mUse2D) {
        double lat2D = par.getValue<double>("MODEL_2D_LATITUDE");
        double lon2D = par.getValue<double>("MODEL_2D_LONGITUDE");
        RDCol3 rtpG;
        rtpG(0) = 1.;
        rtpG(1) = XMath::lat2Theta(lat2D, 0.);
        rtpG(2) = XMath::lon2Phi(lon2D);
        const RDCol3 &rtpS = XMath::rotateGlob2Src(rtpG, srcLat, srcLon, srcDep);
        mPhi2D = rtpS(2);
    } 
    
    // slice plots
    SlicePlot::buildInparam(mSlicePlots, par, this);
    if (mSlicePlots.size() > 0 && XMPI::root()) {
        std::fstream fs;
        // node coordinates
        fs.open(Parameters::sOutputDirectory + "/plots/mesh_coordinates.txt", std::fstream::out);
        for (int i = 0; i < exModel->getNumNodes(); i++) 
            fs << exModel->getNodalS(i) << " " << exModel->getNodalZ(i) << std::endl; 
        fs.close();
        // element connectivity
        const std::vector<std::array<int, 4>> &con = exModel->getConnectivity();
        fs.open(Parameters::sOutputDirectory + "/plots/mesh_connectivity.txt", std::fstream::out);
        for (int i = 0; i < exModel->getNumQuads(); i++) 
            fs << con[i][0] << " " << con[i][1] << " " << con[i][2] << " " << con[i][3] << std::endl;
        fs.close();
    }
}

void Mesh::buildUnweighted() {
    int nElemGlobal = mExModel->getNumQuads(); 
    // balance solid and fluid separately
    DecomposeOption option(nElemGlobal, true, 0., 1, mDDPar->mNPartMetis, false);
    for (int iquad = 0; iquad < nElemGlobal; iquad++) {
        if (mExModel->getElementalVariables().at("fluid")[iquad] > .5) 
            option.mElemWeights2[iquad] = 1.; // fluid
        else 
            option.mElemWeights1[iquad] = 1.; // solid
    }
    XTimer::begin("Build Local", 1);
    buildLocal(option);
    XTimer::end("Build Local", 1);
    // slice plots
    XTimer::begin("Plot at Unweighted Phase", 1);
    for (const auto &sp: mSlicePlots) sp->plotUnweighted();
    XTimer::end("Plot at Unweighted Phase", 1);
}

double Mesh::getDeltaT() const {
    double dt = 1e100;
    for (int i = 0; i < getNumQuads(); i++) 
        dt = std::min(dt, mQuads[i]->getDeltaT());
    return XMPI::min(dt);
}

void Mesh::setAttBuilder(const AttBuilder *attBuild) {
    mAttBuilder = attBuild;
}

void Mesh::buildWeighted() {
    DecomposeOption measured;
    XTimer::begin("Measure", 1);
    measure(measured);
    XTimer::end("Measure", 1);
    
    XTimer::begin("Build Local", 1);
    buildLocal(measured);
    XTimer::end("Build Local", 1);
    // slice plots
    XTimer::begin("Plot at Weighted Phase", 1);
    for (const auto &sp: mSlicePlots) sp->plotWeighted();
    XTimer::end("Plot at Weighted Phase", 1);
}

void Mesh::release(Domain &domain) {
    XTimer::begin("Release Points", 2);
    for (const auto &point: mGLLPoints) point->release(domain);
    XTimer::end("Release Points", 2);
    
    XTimer::begin("Release Elements", 2);
    for (int iloc = 0; iloc < getNumQuads(); iloc++) {
        int etag = mQuads[iloc]->release(domain, mLocalElemToGLL[iloc], mAttBuilder);
        mQuads[iloc]->setElementTag(etag);
    }
    XTimer::end("Release Elements", 2);
    
    // set messaging 
    MessagingBuffer *buf = new MessagingBuffer();
    for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
        int sz_total = 0;
        int npoint = mMsgInfo->mNLocalPoints[i];
        for (int j = 0; j < npoint; j++) {
            int pTag = mMsgInfo->mILocalPoints[i][j];
            sz_total += domain.getPoint(pTag)->sizeComm();
        }
        buf->mBufferSend.push_back(CColX(sz_total));
        buf->mBufferRecv.push_back(CColX(sz_total));
    }    
    MessagingInfo *msg = new MessagingInfo(*mMsgInfo);
    domain.setMessaging(msg, buf);
    
    // set learn parameters
    domain.setLearnParameters(new LearnParameters(*mLearnPar));
}

double Mesh::computeRadiusRef(double depth, double lat, double lon) const {
    // geocentric 
    double theta = XMath::lat2Theta(lat, depth);
    double phi = XMath::lon2Phi(lon);
    double router = mExModel->getROuter();
    
    // 2D mode
    if (mUse2D) {
        RDCol3 rtpG, rtpS;
        rtpG(0) = 1.;
        rtpG(1) = theta;
        rtpG(2) = phi;
        rtpS = XMath::rotateGlob2Src(rtpG, mSrcLat, mSrcLon, mSrcDep);
        // enforced azimuth
        rtpS(2) = mPhi2D;
        rtpG = XMath::rotateSrc2Glob(rtpS, mSrcLat, mSrcLon, mSrcDep);
        theta = rtpG(1);
        phi = rtpG(2);
    }
    
    // surface receivers
    if (depth < tinyDouble) return router;
    
    // target
    double R = computeRPhysical(router, theta, phi) - depth;
    double distTol = std::min((double)tinySingle, mExModel->getDistTolerance() * tinySingle);
    
    // computeRPhysical is monotonically increasing, so we use Divide-and-Conquer method
    // initial guess = router - depth
    double current = router - depth;
    double upper = router;
    double lower = 0.;
    int maxIter = 10000;
    int iter = 0;
    while (iter++ <= maxIter) {
        double diff = computeRPhysical(current, theta, phi) - R;
        if (std::abs(diff) < distTol) return current;
        if (diff > 0.) 
            upper = current;
        else 
            lower = current;
        current = .5 * (lower + upper);
    }
    throw std::runtime_error("Mesh::computeRadiusRef || Failed to find reference radius.");
}

double Mesh::computeRPhysical(double r, double theta, double phi) const {
    double deltaR = 0.;
    for (const auto &g3D: mGeometric3D) 
        deltaR += g3D->getDeltaR(r, theta, phi, r);
    return r + deltaR;
}

void Mesh::buildLocal(const DecomposeOption &option) {
    // destroy existent
    destroy();
    
    // domain decomposition
    XTimer::begin("Domain Decomposition", 2);
    int nGllLocal;
    IColX procMask;
    mMsgInfo = new MessagingInfo();
    Connectivity con_global(mExModel->getConnectivity());
    con_global.decompose(option, nGllLocal, mLocalElemToGLL, *mMsgInfo, procMask);
    XTimer::end("Domain Decomposition", 2);
    
    // form empty points 
    XTimer::begin("Generate Points", 2);
    mGLLPoints.reserve(nGllLocal);
    for (int i = 0; i < nGllLocal; i++) 
        mGLLPoints.push_back(new GLLPoint());
    XTimer::end("Generate Points", 2);
    
    // quads 
    XTimer::begin("Generate Quads", 2);
    double s_max, s_min, z_max, z_min;
    mSMax = mZMax = -1e30;
    mSMin = mZMin = 1e30;
    mQuads.reserve(mLocalElemToGLL.size());
    for (int iquad = 0; iquad < mExModel->getNumQuads(); iquad++) {
        if (procMask(iquad)) {
            // 1D Quad
            Quad *quad = new Quad(*mExModel, iquad, *mNrField);
            // 3D model
            for (int j = 0; j < mVolumetric3D.size(); j++) 
                quad->addVolumetric3D(*(mVolumetric3D[j]), mSrcLat, mSrcLon, mSrcDep, mPhi2D);
            for (int j = 0; j < mGeometric3D.size(); j++) 
                quad->addGeometric3D(*(mGeometric3D[j]), mSrcLat, mSrcLon, mSrcDep, mPhi2D);
            if (mOceanLoad3D != 0) quad->setOceanLoad3D(*mOceanLoad3D, mSrcLat, mSrcLon, mSrcDep, mPhi2D);    
            quad->finishModel3D();
            // spatial range
            quad->getSpatialRange(s_max, s_min, z_max, z_min);
            mSMax = std::max(mSMax, s_max);
            mSMin = std::min(mSMin, s_min);
            mZMax = std::max(mZMax, z_max);
            mZMin = std::min(mZMin, z_min);
            // add to local build
            mQuads.push_back(quad);
        }
    }
    XTimer::end("Generate Quads", 2);
    
    // setup GLL points
    XTimer::begin("Setup Points", 2);
    for (int iloc = 0; iloc < mLocalElemToGLL.size(); iloc++) {
        mQuads[iloc]->setupGLLPoints(mGLLPoints, mLocalElemToGLL[iloc], mExModel->getDistTolerance());
    } 
    XTimer::end("Setup Points", 2);
    
    /////////////////////////////// assemble mass and normal ///////////////////////////////
    XTimer::begin("Assemble Mass", 2);
    // mpi buffer
    std::vector<RDMatXX> bufferGLLSend;
    std::vector<RDMatXX> bufferGLLRecv;
    for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
        int nr_max = -1;
        int npoint = mMsgInfo->mNLocalPoints[i];
        for (int j = 0; j < npoint; j++) {
            int pTag = mMsgInfo->mILocalPoints[i][j];
            int nr = mGLLPoints[pTag]->getNr();
            nr_max = std::max(nr, nr_max);
        }
        // nr_max for solid mass
        // nr_max for fluid mass
        // nr_max * 3 for surface normal
        // nr_max * 3 for solid-fluid normal
        // 1 for reference count
        bufferGLLSend.push_back(RDMatXX::Zero(nr_max * 8 + 1, npoint));
        bufferGLLRecv.push_back(RDMatXX::Zero(nr_max * 8 + 1, npoint));
    }    
    
    // feed buffer
    for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
        for (int j = 0; j < mMsgInfo->mNLocalPoints[i]; j++) {
            int pTag = mMsgInfo->mILocalPoints[i][j];
            mGLLPoints[pTag]->feedBuffer(bufferGLLSend[i], j);
        }
    }
    
    // send and recv
    for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
        mMsgInfo->mReqSend[i] = XMPI::isend(mMsgInfo->mIProcComm[i], bufferGLLSend[i]);
        mMsgInfo->mReqRecv[i] = XMPI::irecv(mMsgInfo->mIProcComm[i], bufferGLLRecv[i]);
    }
    
    // wait recv
    XMPI::wait_all(mMsgInfo->mReqRecv.begin(), mMsgInfo->mReqRecv.end());
    
    // extract buffer 
    for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
        for (int j = 0; j < mMsgInfo->mNLocalPoints[i]; j++) {
            int pTag = mMsgInfo->mILocalPoints[i][j];
            mGLLPoints[pTag]->extractBuffer(bufferGLLRecv[i], j);
        }
    }
    
    // wait send 
    XMPI::wait_all(mMsgInfo->mReqSend.begin(), mMsgInfo->mReqSend.end());
    
    XTimer::end("Assemble Mass", 2);
}

void Mesh::destroy() {
    // points
    for (const auto &point: mGLLPoints) delete point;
    mGLLPoints.clear();
    // quads
    for (const auto &quad: mQuads) delete quad;
    mQuads.clear();
    // l2g mapping
    mLocalElemToGLL.clear();
    // message
    if (mMsgInfo) {
        delete mMsgInfo;
        mMsgInfo = 0;
    }
}

void Mesh::measure(DecomposeOption &measured) {
    // a temp Domain
    Domain domain;
    release(domain);
    
    // user clock resolution
    bool useClockUser = false;
    double clockFactor = useClockUser ? 1e2 : 1e4;
    double clockResolution = XMath::getClockResolution(useClockUser);
    // minimum number of steps to be measured 
    int minStep = 5;
    // how may elements of the same kind will be measured
    int nMeasureSameKind = 3;
    
    ////////// measure elements //////////
    XTimer::begin("Measure Elements", 2);
    // initialize with zero weights
    int nElemGlobal = mExModel->getNumQuads();
    RDColX eWgtEle = RDColX::Zero(nElemGlobal);
    IColX eCommSize = IColX::Zero(nElemGlobal);
    // create library
    std::map<std::string, double> elemCostLibrary;
    for (int iloc = 0; iloc < getNumQuads(); iloc++) {
        int elemTag = mQuads[iloc]->getElementTag();
        Element *elem = domain.getElement(elemTag);
        // get cost signature
        std::string coststr = elem->costSignature();
        // insert to library
        elemCostLibrary.insert(std::pair<std::string, double>(coststr, -1.));
        // perform measurement only if it is new (measure = -1.)
        if (elemCostLibrary.at(coststr) < 0.) {
            // find how may steps are needed to use USER clock
            double wall = elem->measure(minStep, false);
            int nstep = std::max(minStep, (int)(clockResolution * clockFactor / wall) + 1);
            elemCostLibrary.at(coststr) = elem->measure(nstep, useClockUser);
            // find elements with the same signature
            int sameKindFound = 0;
            for (int jloc = iloc + 1; jloc < getNumQuads(); jloc++) {
                int elemTagOther = mQuads[jloc]->getElementTag();
                Element *elemOther = domain.getElement(elemTagOther);
                if (elemOther->costSignature() == coststr) {
                    // use minimum
                    elemCostLibrary.at(coststr) = std::min(elemOther->measure(nstep, useClockUser), 
                        elemCostLibrary.at(coststr));
                    if (++sameKindFound == nMeasureSameKind) break;
                }
            }
        }
    }
    XTimer::end("Measure Elements", 2);
    
    // uniform measurements across procs
    XTimer::begin("Bcast Element Costs", 2);
    std::vector<std::map<std::string, double>> all_elemCostLibrary = XMPI::all_gather(elemCostLibrary);
    std::map<std::string, double> elemCostLibraryGlobal;
    for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
        const std::map<std::string, double> &iproc_lib = all_elemCostLibrary[iproc];
        for (auto it = iproc_lib.begin(); it != iproc_lib.end(); it++) {
            elemCostLibraryGlobal.insert(*it);
            elemCostLibraryGlobal.at(it->first) = std::min(it->second,
                elemCostLibraryGlobal.at(it->first));
        }
    }
    // read library
    for (int iloc = 0; iloc < getNumQuads(); iloc++) {
        int elemTag = mQuads[iloc]->getElementTag();
        Element *elem = domain.getElement(elemTag);
        // read library
        double measure = elemCostLibraryGlobal.at(elem->costSignature());
        // assign to element
        int quadTag = mQuads[iloc]->getQuadTag();
        eWgtEle(quadTag) = measure;
        // commnunication size    
        eCommSize(quadTag) = elem->sizeComm();
    }
    XTimer::end("Bcast Element Costs", 2);
    

    ////////// measure points //////////
    XTimer::begin("Measure Points", 2);
    // initialize with zero weights
    int ngll = mGLLPoints.size();
    RDColX pWgt = RDColX::Zero(ngll);
    // create library
    std::map<std::string, double> pointCostLibrary;
    for (int ip = 0; ip < ngll; ip++) {
        Point *point = domain.getPoint(ip);
        // get cost signature
        std::string coststr = point->costSignature();
        // insert to library
        pointCostLibrary.insert(std::pair<std::string, double>(coststr, -1.));
        // perform measurement only if it is new (measure = -1.)
        if (pointCostLibrary.at(coststr) < 0.) {
            // find how may steps are needed to use USER clock
            double wall = point->measure(minStep, false);
            int nstep = std::max(minStep, (int)(clockResolution * clockFactor / 10. / wall) + 1);
            pointCostLibrary.at(coststr) = point->measure(nstep, useClockUser);
            // find points with the same signature
            int sameKindFound = 0;
            for (int jp = ip + 1; jp < ngll; jp++) {
                Point *pointOther = domain.getPoint(jp);
                if (pointOther->costSignature() == coststr) {
                    pointCostLibrary.at(coststr) = std::min(pointOther->measure(nstep, useClockUser),
                        pointCostLibrary.at(coststr));
                    if (++sameKindFound == nMeasureSameKind) break;
                }
            }
        }
    }
    XTimer::end("Measure Points", 2);
    
    // uniform measurements across procs
    XTimer::begin("Bcast Point Costs", 2);
    std::vector<std::map<std::string, double>> all_pointCostLibrary = XMPI::all_gather(pointCostLibrary);
    std::map<std::string, double> pointCostLibraryGlobal;
    for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
        const std::map<std::string, double> &iproc_lib = all_pointCostLibrary[iproc];
        for (auto it = iproc_lib.begin(); it != iproc_lib.end(); it++) {
            pointCostLibraryGlobal.insert(*it);
            pointCostLibraryGlobal.at(it->first) = std::min(it->second,
                pointCostLibraryGlobal.at(it->first));
        }
    }
    // read library
    for (int ip = 0; ip < ngll; ip++) {
        Point *point = domain.getPoint(ip);
        // read library
        double measure = pointCostLibraryGlobal.at(point->costSignature());
        // assign to point
        pWgt(ip) = measure / mGLLPoints[ip]->getReferenceCount();
    }
    XTimer::end("Bcast Point Costs", 2);
    
    XTimer::begin("Finish Measurement", 2);    
    // add point weights to elements
    RDColX eWgtPnt = RDColX::Zero(nElemGlobal);
    for (int iloc = 0; iloc < getNumQuads(); iloc++) {
        int quadTag = mQuads[iloc]->getQuadTag();
        for (int ipol = 0; ipol <= nPol; ipol++) 
            for (int jpol = 0; jpol <= nPol; jpol++) 
                eWgtPnt(quadTag) += pWgt(mLocalElemToGLL[iloc](ipol, jpol));
    }
    
    // sum up
    XMPI::sumEigen(eWgtEle);
    XMPI::sumEigen(eWgtPnt);
    XMPI::sumEigen(eCommSize);
    
    // create option 
    if (mDDPar->mBalanceEP) {
        measured = DecomposeOption(nElemGlobal, true, 0., 0, mDDPar->mNPartMetis, mDDPar->mCommVolMetis);
        for (int i = 0; i < nElemGlobal; i++) {
            measured.mElemWeights1[i] = eWgtEle(i);
            measured.mElemWeights2[i] = eWgtPnt(i);
            measured.mElemCommSize[i] = eCommSize(i);
        }
    } else {
        measured = DecomposeOption(nElemGlobal, false, 0., 0, mDDPar->mNPartMetis, mDDPar->mCommVolMetis);
        for (int i = 0; i < nElemGlobal; i++) {
            measured.mElemWeights1[i] = eWgtEle(i) + eWgtPnt(i);
            measured.mElemCommSize[i] = eCommSize(i);
        }
    }
    
    // report
    if (XMPI::root() && mDDPar->mReportMeasure) {
        std::string fname = Parameters::sOutputDirectory + "/develop/measured_costs.txt";
        std::fstream fs(fname, std::fstream::out);
        fs << "*** Element Types ***" << std::endl;
        for (auto it = elemCostLibraryGlobal.begin(); it != elemCostLibraryGlobal.end(); it++) 
            fs << it->first << "    " << it->second << std::endl;
        fs << std::endl << "*** Point Types ***" << std::endl;    
        for (auto it = pointCostLibraryGlobal.begin(); it != pointCostLibraryGlobal.end(); it++) 
            fs << it->first << "    " << it->second << std::endl;    
        fs.close();    
    }
    
    // plot 
    XTimer::begin("Plot during Cost Measurements", 2);
    for (const auto &sp: mSlicePlots) {
        sp->plotEleType(domain);
        sp->plotMeasured(eWgtEle + eWgtPnt);
    }
    XTimer::end("Plot during Cost Measurements", 2);
    
    XTimer::end("Finish Measurement", 2);    
}

void Mesh::test() {
    // a temp Domain
    Domain domain;
    release(domain);
    domain.test();
}

Mesh::DDParameters::DDParameters(const Parameters &par) {
    mBalanceEP = par.getValue<bool>("DD_BALANCE_ELEMENT_POINT");
    mNPartMetis = par.getValue<int>("DD_NPART_METIS");
    mCommVolMetis = par.getValue<bool>("DD_COMM_VOL_METIS");
    mReportMeasure = par.getValue<bool>("DEVELOP_MEASURED_COSTS");
    if (mNPartMetis <= 0) mNPartMetis = 10;
}

