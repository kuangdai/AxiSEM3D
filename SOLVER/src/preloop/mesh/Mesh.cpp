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
#include "Geodesy.h"
#include "Geometric3D.h"
#include "NuWisdom.h"
#include "NrField.h"

#include "MultilevelTimer.h"
#include "SlicePlot.h"
#include <fstream>
#include <cfloat>

Mesh::~Mesh() {
    destroy(); // local build
    delete mDDPar;
    delete mLearnPar;
    for (const auto &sp: mSlicePlots) {
        delete sp;
    }    
}

Mesh::Mesh(const ExodusModel *exModel, const NrField *nrf, 
    double srcLat, double srcLon, double srcDep, const Parameters &par, int verbose):
mExModel(exModel), mNrField(nrf), mSrcLat(srcLat), mSrcLon(srcLon), mSrcDep(srcDep) {
    mAttBuilder = 0;
    mMsgInfo = 0;
    mOceanLoad3D = 0;
    mDDPar = new DDParameters(par);
    mLearnPar = new LearnParameters(par);
    
    // 2D mode
    std::string mode2d = par.getValue<std::string>("MODEL_2D_MODE");
    if (boost::iequals(mode2d, "off")) {
        mPhi2D = -DBL_MAX;
    } else if (boost::iequals(mode2d, "geographic")) {
        double lat2D = par.getValue<double>("MODEL_2D_LATITUDE");
        double lon2D = par.getValue<double>("MODEL_2D_LONGITUDE");
        RDCol3 rtpG;
        rtpG(0) = 1.;
        rtpG(1) = Geodesy::lat2Theta_d(lat2D, 0.);
        rtpG(2) = Geodesy::lon2Phi(lon2D);
        const RDCol3 &rtpS = Geodesy::rotateGlob2Src(rtpG, srcLat, srcLon, srcDep);
        mPhi2D = rtpS(2);
    } else if (boost::iequals(mode2d, "source-centered")) {
         mPhi2D = par.getValue<double>("MODEL_2D_AZIMUTH") * degree;
         while (mPhi2D < 0.) {
             mPhi2D += 2. * pi;
         }
         while (mPhi2D > 2. * pi) {
             mPhi2D -= 2. * pi;
         }
    } else {
         throw std::runtime_error("Mesh::Mesh || Invalid input for MODEL_2D_MODE.");
    }
    
    // slice plots
    SlicePlot::buildInparam(mSlicePlots, par, this, verbose);
    if (mSlicePlots.size() > 0 && XMPI::root()) {
        std::fstream fs;
        // node coordinates
        fs.open(Parameters::sOutputDirectory + "/plots/mesh_coordinates.txt", std::fstream::out);
        for (int i = 0; i < exModel->getNumNodes(); i++) {
            fs << exModel->getNodalS(i) << " " << exModel->getNodalZ(i) << std::endl; 
        }
        fs.close();
        // element connectivity
        fs.open(Parameters::sOutputDirectory + "/plots/mesh_connectivity.txt", std::fstream::out);
        fs << exModel->getConnectivity() << std::endl;
        fs.close();
    }
}

void Mesh::buildUnweighted() {
    // balance by Nr
    DecomposeOption option;
    option.mProcInterval = mDDPar->mProcInterval;
    option.mNCutsPerProc = mDDPar->mNCutsPerProc;
    int nElemGlobal = mExModel->getNumQuads(); 
    option.mElemWeights = RDColX::Zero(nElemGlobal);
    for (int iquad = 0; iquad < nElemGlobal; iquad++) {
        int inode = mExModel->getConnectivity()(iquad, 0);
        RDCol2 sz;
        sz(0) = mExModel->getNodalS(inode);
        sz(1) = mExModel->getNodalZ(inode);
        option.mElemWeights(iquad) = mNrField->getNrAtPoint(sz) * 1.;
    }
    MultilevelTimer::begin("Build Local", 1);
    buildLocal(option);
    MultilevelTimer::end("Build Local", 1);
    // slice plots
    MultilevelTimer::begin("Plot at Unweighted Phase", 1);
    for (const auto &sp: mSlicePlots) sp->plotUnweighted();
    MultilevelTimer::end("Plot at Unweighted Phase", 1);
}

double Mesh::getDeltaT() const {
    double dt = DBL_MAX;
    for (int i = 0; i < getNumQuads(); i++) {
        dt = std::min(dt, mQuads[i]->getDeltaT());
    }
    return XMPI::min(dt);
}

void Mesh::setAttBuilder(const AttBuilder *attBuild) {
    mAttBuilder = attBuild;
}

void Mesh::buildWeighted() {
    DecomposeOption measured;
    measured.mProcInterval = mDDPar->mProcInterval;
    measured.mNCutsPerProc = mDDPar->mNCutsPerProc;
    MultilevelTimer::begin("Measure", 1);
    measure(measured);
    MultilevelTimer::end("Measure", 1);
    
    MultilevelTimer::begin("Build Local", 1);
    buildLocal(measured);
    MultilevelTimer::end("Build Local", 1);
    // slice plots
    MultilevelTimer::begin("Plot at Weighted Phase", 1);
    for (const auto &sp: mSlicePlots) {
        sp->plotWeighted();
    }
    MultilevelTimer::end("Plot at Weighted Phase", 1);
}

void Mesh::release(Domain &domain) {
    MultilevelTimer::begin("Release Points", 2);
    for (const auto &point: mGLLPoints) {
        point->release(domain);
    }
    MultilevelTimer::end("Release Points", 2);
    
    MultilevelTimer::begin("Release Elements", 2);
    for (int iloc = 0; iloc < getNumQuads(); iloc++) {
        int etag = mQuads[iloc]->release(domain, mLocalElemToGLL[iloc], mAttBuilder);
        mQuads[iloc]->setElementTag(etag);
    }
    MultilevelTimer::end("Release Elements", 2);
    
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
    double theta = Geodesy::lat2Theta_d(lat, depth);
    double phi = Geodesy::lon2Phi(lon);
    double router = mExModel->getROuter();
    
    // 2D mode
    if (mPhi2D < -DBL_MAX * .9) {
        RDCol3 rtpG, rtpS;
        rtpG(0) = 1.;
        rtpG(1) = theta;
        rtpG(2) = phi;
        rtpS = Geodesy::rotateGlob2Src(rtpG, mSrcLat, mSrcLon, mSrcDep);
        // enforced azimuth
        rtpS(2) = mPhi2D;
        rtpG = Geodesy::rotateSrc2Glob(rtpS, mSrcLat, mSrcLon, mSrcDep);
        theta = rtpG(1);
        phi = rtpG(2);
    }
    
    // surface receivers
    if (depth < tinyDouble) {
        return router;
    }
    
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
        if (std::abs(diff) < distTol) {
            return current;
        }
        if (diff > 0.) {
            upper = current;
        } else {
            lower = current;
        }
        current = .5 * (lower + upper);
    }
    throw std::runtime_error("Mesh::computeRadiusRef || Failed to find reference radius.");
}

double Mesh::computeRPhysical(double r, double theta, double phi) const {
    double deltaR = 0.;
    for (const auto &g3D: mGeometric3D) {
        deltaR += g3D->getDeltaR(r, theta, phi, r);
    }
    return r + deltaR;
}

// # include "NetCDF_Writer.h"
void Mesh::buildLocal(const DecomposeOption &option) {
    // destroy existent
    destroy();
    
    // domain decomposition
    MultilevelTimer::begin("Domain Decomposition", 2);
    int nGllLocal;
    IColX procMask;
    mMsgInfo = new MessagingInfo();
    Connectivity con_global(mExModel->getConnectivity());
    con_global.decompose(option, nGllLocal, mLocalElemToGLL, *mMsgInfo, procMask);
    MultilevelTimer::end("Domain Decomposition", 2);
    
    // form empty points 
    MultilevelTimer::begin("Generate Points", 2);
    mGLLPoints.reserve(nGllLocal);
    for (int i = 0; i < nGllLocal; i++) {
        mGLLPoints.push_back(new GLLPoint());
    }
    MultilevelTimer::end("Generate Points", 2);
    
    // quads 
    MultilevelTimer::begin("Generate Quads", 2);
    double s_max, s_min, z_max, z_min;
    mSMax = mZMax = -1e30;
    mSMin = mZMin = 1e30;
    mQuads.reserve(mLocalElemToGLL.size());
    for (int iquad = 0; iquad < mExModel->getNumQuads(); iquad++) {
        if (procMask(iquad)) {
            // 1D Quad
            Quad *quad = new Quad(*mExModel, iquad, *mNrField);
            // 3D model
            quad->addVolumetric3D(mVolumetric3D, mSrcLat, mSrcLon, mSrcDep, mPhi2D);
            quad->addGeometric3D(mGeometric3D, mSrcLat, mSrcLon, mSrcDep, mPhi2D);
            if (mOceanLoad3D != 0) {
                quad->setOceanLoad3D(*mOceanLoad3D, mSrcLat, mSrcLon, mSrcDep, mPhi2D);
            }    
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
    MultilevelTimer::end("Generate Quads", 2);
    
    // setup GLL points
    MultilevelTimer::begin("Setup Points", 2);
    for (int iloc = 0; iloc < mLocalElemToGLL.size(); iloc++) {
        mQuads[iloc]->setupGLLPoints(mGLLPoints, mLocalElemToGLL[iloc], mExModel->getDistTolerance());
    } 
    MultilevelTimer::end("Setup Points", 2);
    
    // RDMatXX_RM sz = RDMatXX::Zero(nGllLocal, 2);
    // for (int i = 0; i < nGllLocal; i++) {
    //     sz.row(i) = mGLLPoints[i]->getCoords().transpose();
    // }
    // NetCDF_Writer nc;
    // nc.open(Parameters::sOutputDirectory + "/plots/gllsz.nc", true);
    // std::vector<size_t> dims;
    // dims.push_back(nGllLocal);
    // dims.push_back(2);
    // nc.defineVariable<double>("gllsz", dims);
    // nc.writeVariableWhole("gllsz", sz);
    // exit(0);
    
    /////////////////////////////// assemble mass and normal ///////////////////////////////
    MultilevelTimer::begin("Assemble Mass", 2);
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
        XMPI::isendDouble(mMsgInfo->mIProcComm[i], bufferGLLSend[i], mMsgInfo->mReqSend[i]);
        XMPI::irecvDouble(mMsgInfo->mIProcComm[i], bufferGLLRecv[i], mMsgInfo->mReqRecv[i]);
    }
    
    // wait recv
    XMPI::wait_all(mMsgInfo->mReqRecv.size(), mMsgInfo->mReqRecv.data());
    
    // extract buffer 
    for (int i = 0; i < mMsgInfo->mNProcComm; i++) {
        for (int j = 0; j < mMsgInfo->mNLocalPoints[i]; j++) {
            int pTag = mMsgInfo->mILocalPoints[i][j];
            mGLLPoints[pTag]->extractBuffer(bufferGLLRecv[i], j);
        }
    }
    
    // wait send 
    XMPI::wait_all(mMsgInfo->mReqSend.size(), mMsgInfo->mReqSend.data());
    
    MultilevelTimer::end("Assemble Mass", 2);
}

void Mesh::destroy() {
    // points
    for (const auto &point: mGLLPoints) {
        delete point;
    }
    mGLLPoints.clear();
    // quads
    for (const auto &quad: mQuads) {
        delete quad;
    }
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
    double clockFactor = 1e4;
    double clockResolution = MyBoostTimer::getClockResolution();
    // minimum number of steps to be measured 
    int minStep = 5;
    // how may elements of the same kind will be measured
    int nMeasureSameKind = 3;
    
    ////////// measure elements //////////
    MultilevelTimer::begin("Measure Elements", 2);
    // initialize with zero weights
    int nElemGlobal = mExModel->getNumQuads();
    RDColX eWgtEle = RDColX::Zero(nElemGlobal);
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
            double wall = elem->measure(minStep);
            int nstep = std::max(minStep, (int)(clockResolution * clockFactor / wall) + 1);
            elemCostLibrary.at(coststr) = elem->measure(nstep);
            // find elements with the same signature
            int sameKindFound = 0;
            for (int jloc = iloc + 1; jloc < getNumQuads(); jloc++) {
                int elemTagOther = mQuads[jloc]->getElementTag();
                Element *elemOther = domain.getElement(elemTagOther);
                if (elemOther->costSignature() == coststr) {
                    // use minimum
                    elemCostLibrary.at(coststr) = std::min(elemOther->measure(nstep), 
                        elemCostLibrary.at(coststr));
                    if (++sameKindFound == nMeasureSameKind) {
                        break;
                    }
                }
            }
        }
    }
    MultilevelTimer::end("Measure Elements", 2);
    
    // uniform measurements across procs
    MultilevelTimer::begin("Bcast Element Costs", 2);
    std::vector<std::map<std::string, double>> all_elemCostLibrary;
    XMPI::gather(elemCostLibrary, all_elemCostLibrary, MPI_DOUBLE, true);
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
    }
    // sum up
    XMPI::sumEigenDouble(eWgtEle);
    MultilevelTimer::end("Bcast Element Costs", 2);
    

    ////////// measure points //////////
    MultilevelTimer::begin("Measure Points", 2);
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
            double wall = point->measure(minStep);
            int nstep = std::max(minStep, (int)(clockResolution * clockFactor / 10. / wall) + 1);
            pointCostLibrary.at(coststr) = point->measure(nstep);
            // find points with the same signature
            int sameKindFound = 0;
            for (int jp = ip + 1; jp < ngll; jp++) {
                Point *pointOther = domain.getPoint(jp);
                if (pointOther->costSignature() == coststr) {
                    pointCostLibrary.at(coststr) = std::min(pointOther->measure(nstep),
                        pointCostLibrary.at(coststr));
                    if (++sameKindFound == nMeasureSameKind) {
                        break;
                    }
                }
            }
        }
    }
    MultilevelTimer::end("Measure Points", 2);
    
    // uniform measurements across procs
    MultilevelTimer::begin("Bcast Point Costs", 2);
    std::vector<std::map<std::string, double>> all_pointCostLibrary;
    XMPI::gather(pointCostLibrary, all_pointCostLibrary, MPI_DOUBLE, true);
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
        pWgt(ip) = pointCostLibraryGlobal.at(point->costSignature());
    }
    MultilevelTimer::end("Bcast Point Costs", 2);
    
    // report
    MultilevelTimer::begin("Report Measurements", 2);
    if (XMPI::root() && mDDPar->mReportMeasure) {
        std::string fname = Parameters::sOutputDirectory + "/develop/measured_costs.txt";
        std::fstream fs(fname, std::fstream::out);
        fs << "*** Element Types ***" << std::endl;
        for (auto it = elemCostLibraryGlobal.begin(); it != elemCostLibraryGlobal.end(); it++) {
            fs << it->first << "    " << it->second << std::endl;
        }
        fs << std::endl << "*** Point Types ***" << std::endl;    
        for (auto it = pointCostLibraryGlobal.begin(); it != pointCostLibraryGlobal.end(); it++) {
            fs << it->first << "    " << it->second << std::endl;    
        }
        fs.close();    
    }
    MultilevelTimer::end("Report Measurements", 2);
    
    // create option 
    MultilevelTimer::begin("Create Weights", 2);    
    // recast point weights element-wise, considering reference count on edges
    RDColX eWgtPnt = RDColX::Zero(nElemGlobal);
    for (int iloc = 0; iloc < getNumQuads(); iloc++) {
        int quadTag = mQuads[iloc]->getQuadTag();
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                int ip = mLocalElemToGLL[iloc](ipol, jpol);
                eWgtPnt(quadTag) += pWgt(ip) / mGLLPoints[ip]->getReferenceCount();
            }
        }
    }
    XMPI::sumEigenDouble(eWgtPnt);
    // sum up element and point weights
    measured.mElemWeights = eWgtEle + eWgtPnt;
    MultilevelTimer::end("Create Weights", 2);    
    
    // plot 
    MultilevelTimer::begin("Plot during Cost Measurements", 2);
    for (const auto &sp: mSlicePlots) {
        sp->plotEleType(domain);
        sp->plotMeasured(measured.mElemWeights);
    }
    MultilevelTimer::end("Plot during Cost Measurements", 2);    
}

void Mesh::test() {
    // a temp Domain
    Domain domain;
    release(domain);
    domain.test();
}

Mesh::DDParameters::DDParameters(const Parameters &par) {
    mReportMeasure = par.getValue<bool>("DEVELOP_MEASURED_COSTS");
    mProcInterval = par.getValue<int>("DD_PROC_INTERVAL");
    mNCutsPerProc = par.getValue<int>("DD_NCUTS_PER_PROC");
    if (mProcInterval <= 0) {
        mProcInterval = 1;
    }
    if (mNCutsPerProc <= 0) {
        mNCutsPerProc = 1;
    }
}

int Mesh::getMaxNr() const {
    int maxNr = -1;
    for (int i = 0; i < getNumQuads(); i++) {
        maxNr = std::max(maxNr, mQuads[i]->getNr());
    }
    return XMPI::max(maxNr);
}
