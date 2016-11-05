// axisem.cpp
// created by Kuangdai on 26-Mar-2016 
// main of AxiSEM3D

#include "axisem.h"
#include "XTimer.h"

int axisem_main(int argc, char *argv[]) {
    
    try {
        
        // variable sets
        PreloopVariables pl;
        SolverVariables sv;
        
        // initialize mpi
        XMPI::initialize(argc, argv);
        
        //////// spectral-element constants
        SpectralConstants::initialize(nPol);  
        
        //////// input parameters 
        int verbose;
        Parameters::buildInparam(pl.mParameters, verbose);
        
        // preloop timer
        XTimer::initialize(Parameters::sOutputDirectory + "/preloop_timer.txt", 4);
        if (pl.mParameters->getValue<bool>("DEVELOP_DIAGNOSE_PRELOOP")) XTimer::enable();
        
        //////// exodus model and attenuation parameters 
        XTimer::begin("Exodus", 0);
        ExodusModel::buildInparam(pl.mExodusModel, *(pl.mParameters), pl.mAttParameters, verbose);
        XTimer::end("Exodus", 0);
        
        //////// fourier field 
        XTimer::begin("NrField", 0);
        NrField::buildInparam(pl.mNrField, *(pl.mParameters), pl.mExodusModel->getROuter(), verbose);
        XTimer::end("NrField", 0);
        
        //////// static variables in solver, mainly FFTW
        XTimer::begin("FFTW", 0);
        initializeSolverStatic(*(pl.mNrField)); 
        XTimer::end("FFTW", 0);
        
        //////// 3D models 
        XTimer::begin("3D Models", 0);
        Volumetric3D::buildInparam(pl.mVolumetric3D, *(pl.mParameters), verbose);
        Geometric3D::buildInparam(pl.mGeometric3D, *(pl.mParameters), verbose);
        OceanLoad3D::buildInparam(pl.mOceanLoad3D, *(pl.mParameters), verbose);
        XTimer::end("3D Models", 0);
        
        //////// source
        XTimer::begin("Source", 0);
        Source::buildInparam(pl.mSource, *(pl.mParameters), verbose);
        double srcLat = pl.mSource->getLatitude();
        double srcLon = pl.mSource->getLongitude();
        double srcDep = pl.mSource->getDepth();
        XTimer::end("Source", 0);
        
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
        
        //////// dt
        XTimer::begin("DT", 0);
        double dt = pl.mParameters->getValue<double>("TIME_DELTA_T");
        if (dt < tinyDouble) dt = pl.mMesh->getDeltaT();
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

extern void initializeSolverStatic(const NrField &nrf) {
    // fftw
    SolverFFTW::importWisdom();
    int maxNr = ceil(nrf.getMaxNr() * 1.1);
    SolverFFTW_1::initialize(maxNr);
    SolverFFTW_3::initialize(maxNr); 
    SolverFFTW_N3::initialize(maxNr);
    SolverFFTW_N6::initialize(maxNr);
    SolverFFTW_N9::initialize(maxNr);
    SolverFFTW::exportWisdom();
    PreloopFFTW::initialize(maxNr);
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
