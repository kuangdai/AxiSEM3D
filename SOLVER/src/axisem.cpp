// axisem.cpp
// created by Kuangdai on 26-Mar-2016 
// main of AxiSEM3D

#include "axisem.h"
#include "XTimer.h"

#include "XMPI.h"
#include "eigenc.h"
#include "eigenp.h"

int axisem_main(int argc, char *argv[]) {
    
    try {
        
        // variable sets
        PreloopVariables pl;
        SolverVariables sv;
        
        // initialize mpi
        XMPI::initialize(argc, argv);
        
        // std::map<std::string, double> a;
        // if (XMPI::root()) 
        // {
        //     a.insert(std::pair<std::string, double>("a", 1.));
        //     a.insert(std::pair<std::string, double>("b", 2.));
        // }
        // else 
        // {
        //     a.insert(std::pair<std::string, double>("c", 3.));
        //     a.insert(std::pair<std::string, double>("d", 4.));
        // }
        // 
        // std::vector<std::map<std::string, double>> aa;
        // XMPI::gather(a, aa, MPI_DOUBLE, true);
        // 
        // if (XMPI::rank() == 2) {
        //     for (int i = 0; i < aa.size();i++){
        //         for (auto it = aa[i].begin(); it != aa[i].end(); it++) {
        //             std::cout << it->first << " " << it->second<<std::endl;
        //         }
        //         std::cout<<"-------"<< std::endl;
        //     }
        // }
        // 
        // exit(0);
        
        // std::vector<std::string> str;
        // if (XMPI::root()) {
        //     str.push_back("foo__");
        //     str.push_back(" ba");
        // }
        // XMPI::sendRecvVector(0, 1, str, str);
        // if (!XMPI::root()) {
        //     std::cout << str.size() << std::endl;
        //     std::cout << str[0] << std::endl;
        //     std::cout << str[1] << std::endl;
        //     exit(0);
        // }
        
        //////// spectral-element constants
        SpectralConstants::initialize(nPol);  
        
        //////// input parameters 
        int verbose;
        Parameters::buildInparam(pl.mParameters, verbose);
        
        //////// preloop timer
        XTimer::initialize(Parameters::sOutputDirectory + "/develop/preloop_timer.txt", 4);
        if (pl.mParameters->getValue<bool>("DEVELOP_DIAGNOSE_PRELOOP")) XTimer::enable();
        
        //////// exodus model and attenuation parameters 
        XTimer::begin("Exodus", 0);
        ExodusModel::buildInparam(pl.mExodusModel, *(pl.mParameters), pl.mAttParameters, verbose);
        XTimer::end("Exodus", 0);
        
        //////// fourier field 
        XTimer::begin("NrField", 0);
        NrField::buildInparam(pl.mNrField, *(pl.mParameters), pl.mExodusModel->getROuter(), verbose);
        XTimer::end("NrField", 0);
        
        //////// source
        XTimer::begin("Source", 0);
        Source::buildInparam(pl.mSource, *(pl.mParameters), verbose);
        double srcLat = pl.mSource->getLatitude();
        double srcLon = pl.mSource->getLongitude();
        double srcDep = pl.mSource->getDepth();
        XTimer::end("Source", 0);
        
        //////// 3D models 
        XTimer::begin("3D Models", 0);
        Volumetric3D::buildInparam(pl.mVolumetric3D, *(pl.mParameters), pl.mExodusModel, 
            srcLat, srcLon, srcDep, verbose);
        Geometric3D::buildInparam(pl.mGeometric3D, *(pl.mParameters), verbose);
        OceanLoad3D::buildInparam(pl.mOceanLoad3D, *(pl.mParameters), verbose);
        XTimer::end("3D Models", 0);
        
        //////// mesh, phase 1
        // define mesh
        XTimer::begin("Mesh Definition", 0);
        pl.mMesh = new Mesh(pl.mExodusModel, pl.mNrField, srcLat, srcLon, srcDep, *(pl.mParameters));
        pl.mMesh->setVolumetric3D(pl.mVolumetric3D);
        pl.mMesh->setGeometric3D(pl.mGeometric3D);
        pl.mMesh->setOceanLoad3D(pl.mOceanLoad3D);
        XTimer::end("Mesh Definition", 0);
        
        // build unweighted local mesh 
        XTimer::begin("Unweighted Mesh", 0);
        pl.mMesh->buildUnweighted();
        XTimer::end("Unweighted Mesh", 0);
        
        //////// static variables in solver, mainly FFTW
        XTimer::begin("Initialize FFTW", 0);
        initializeSolverStatic(pl.mMesh->getMaxNr()); 
        XTimer::end("Initialize FFTW", 0);
        
        //////// dt
        XTimer::begin("DT", 0);
        double dt = pl.mParameters->getValue<double>("TIME_DELTA_T");
        if (dt < tinyDouble) dt = pl.mMesh->getDeltaT();
        double dt_fact = pl.mParameters->getValue<double>("TIME_DELTA_T_FACTOR");
        if (dt_fact < tinyDouble) dt_fact = 1.0;
        dt *= dt_fact;
        XTimer::end("DT", 0);
        
        //////// attenuation
        XTimer::begin("Attenuation", 0);
        AttBuilder::buildInparam(pl.mAttBuilder, *(pl.mParameters), *(pl.mAttParameters), dt, verbose);
        XTimer::end("Attenuation", 0);
        
        //////// mesh, phase 2
        XTimer::begin("Weighted Mesh", 0);
        pl.mMesh->setAttBuilder(pl.mAttBuilder);
        pl.mMesh->buildWeighted();
        XTimer::end("Weighted Mesh", 0);
        
        //////// mesh test 
        // test positive-definiteness and self-adjointness of stiffness and mass matrices
        // better to turn with USE_DOUBLE 
        // pl.mMesh->test();
        // XMPI::barrier();
        // exit(0);
        
        //////// source time function 
        XTimer::begin("Source Time Function", 0);
        STF::buildInparam(pl.mSTF, *(pl.mParameters), dt, verbose);
        XTimer::end("Source Time Function", 0);
        
        //////// receivers
        XTimer::begin("Receivers", 0);
        ReceiverCollection::buildInparam(pl.mReceivers, 
            *(pl.mParameters), srcLat, srcLon, srcDep, verbose);
        XTimer::end("Receivers", 0);    
        
        //////// computational domain
        XTimer::begin("Computationalion Domain", 0);
        sv.mDomain = new Domain();
        
        // release mesh
        XTimer::begin("Release Mesh", 1);
        pl.mMesh->release(*(sv.mDomain));
        XTimer::end("Release Mesh", 1);
        
        // release source 
        XTimer::begin("Release Source", 1);
        pl.mSource->release(*(sv.mDomain), *(pl.mMesh));
        XTimer::end("Release Source", 1);
        
        // release stf 
        XTimer::begin("Release STF", 1);
        pl.mSTF->release(*(sv.mDomain));
        XTimer::end("Release STF", 1);
        
        // release receivers
        XTimer::begin("Release Receivers", 1);
        pl.mReceivers->release(*(sv.mDomain), *(pl.mMesh));
        XTimer::end("Release Receivers", 1);
        
        // verbose domain 
        XTimer::begin("Verbose", 1);
        if (verbose) XMPI::cout << sv.mDomain->verbose();
        XTimer::end("Verbose", 1);
        XTimer::end("Computationalion Domain", 0);
        
        XTimer::finalize();
        
        //////////////////////// PREPROCESS DONE ////////////////////////
        
        //////// Newmark
        int infoInt = pl.mParameters->getValue<int>("OPTION_LOOP_INFO_INTERVAL");
        int stabInt = pl.mParameters->getValue<int>("OPTION_STABILITY_INTERVAL");
        sv.mNewmark = new Newmark(sv.mDomain, infoInt, stabInt);
        
        //////// final preparations
        // finalize preloop variables before time loop starts
        pl.finalize();
        // forbid matrix allocation in time loop
        #ifndef NDEBUG
            Eigen::internal::set_is_malloc_allowed(false);
        #endif
            
        //////// GoGoGo
        XMPI::barrier();
        sv.mNewmark->solve();
        
        //////// finalize solver
        // solver 
        sv.finalize();
        // static variables in solver
        finalizeSolverStatic();
        
        // finalize mpi 
        XMPI::finalize();
        
    } catch (const std::exception &e) {
        // print exception
        XMPI::cout.setp(XMPI::rank());
        XMPI::printException(e);
        
        // abort program
        // TODO 
        // MPI_Abort is necessary here. Otherwise, if an exception
        // is thrown from one of the procs, deadlock will occur.
        // But the problem is, how we free memories on other procs?!
        XMPI::abort();
    }
    
    return 0;
}

#include "SolverFFTW.h"
#include "SolverFFTW_1.h"
#include "SolverFFTW_3.h"
#include "SolverFFTW_N3.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"
#include "PreloopFFTW.h"
#include "SolidElement.h"
#include "FluidElement.h"

extern void initializeSolverStatic(int maxNr) {
    // fftw
    SolverFFTW::importWisdom();
    SolverFFTW_1::initialize(maxNr);
    SolverFFTW_3::initialize(maxNr); 
    SolverFFTW_N3::initialize(maxNr);
    SolverFFTW_N6::initialize(maxNr);
    SolverFFTW_N9::initialize(maxNr);
    SolverFFTW::exportWisdom();
    // PreloopFFTW::initialize(maxNr);
    // element
    SolidElement::initWorkspace(maxNr / 2);
    FluidElement::initWorkspace(maxNr / 2);
};

extern void finalizeSolverStatic() {
    // fftw
    SolverFFTW_1::finalize();
    SolverFFTW_3::finalize(); 
    SolverFFTW_N3::finalize();
    SolverFFTW_N6::finalize();
    SolverFFTW_N9::finalize();
    PreloopFFTW::finalize();
};
