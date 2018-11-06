// Newmark.cpp
// created by Kuangdai on 8-Apr-2016 
// Newmark time scheme

#include "Newmark.h"
#include "Domain.h"
#include "SourceTimeFunction.h"
#include <sstream>
#include "XMPI.h"
#include "MultilevelTimer.h"

Newmark::Newmark(Domain *&domain, int reportInterval, int checkStabInterval, bool randomDispl):
mDomain(domain), mReportInterval(reportInterval), 
mCheckStabInterval(checkStabInterval), mRandomDispl(randomDispl) {
    if (mReportInterval <= 0) mReportInterval = 100;
    if (mCheckStabInterval <= 0) mCheckStabInterval = mReportInterval;
}

void Newmark::solve(int verbose) const {
    if (verbose) {
        XMPI::cout << XMPI::endl;
        XMPI::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" << XMPI::endl;
        XMPI::cout << "TTTTTTTTTT  NEWMARK TIME LOOP STARTS  TTTTTTTTTT" << XMPI::endl;
        XMPI::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" << XMPI::endl << XMPI::endl;
    }
    
    double t = 0. - mDomain->getSTF().getShift();
    double dt = mDomain->getSTF().getDeltaT();
    int maxStep = mDomain->getSTF().getSize();
    mDomain->resetZero();
    if (mRandomDispl) {
        mDomain->initDisplTinyRandom();
    }
    const double sec2h = 1. / 3600.;
    double elapsed_last = 0.;
    MyBoostTimer timer;
    timer.start();
    
    ////////////////////////// loop //////////////////////////
    for (int tstep = 1; tstep <= maxStep; tstep++) {
        // update to next step
        mDomain->updateNewmark(dt);
        
        // source
        mDomain->applySource(tstep - 1);
        
        // element stiffness
        mDomain->computeStiff();
        
        // solid-fluid coupling
        mDomain->coupleSolidFluid();
        
        // assemble phase 1: feed + send + recv 
        mDomain->assembleStiff(-1);
        
        // record seismograms
        mDomain->record(tstep - 1, t);    
        t += dt;
        
        // check stability
        if (tstep % mCheckStabInterval == 0) {
            mDomain->checkStability(dt, tstep, t);
        }
        
        // screen info    
        if (tstep % mReportInterval == 0 && verbose) {
            double elapsed = timer.elapsed() * sec2h;
            double speed = (elapsed - elapsed_last) / mReportInterval;
            double total = speed * maxStep;
            double left = speed * (maxStep - tstep);
            elapsed_last = elapsed;
            std::stringstream ss;
            int percent = (int)(100. * tstep / maxStep);
            ss << "  SIMULATION TIME / sec     =   " << t << XMPI::endl;
            ss << "  TIME STEP / TOTAL STEPS   =   " << tstep << " / " <<  maxStep << " (" << percent << "%)" << XMPI::endl;
            ss << "  WALLTIME ELAPSED  / h     =   " << elapsed << XMPI::endl; 
            ss << "  WALLTIME REMAINED / h     =   " << left << XMPI::endl; 
            ss << "  WALLTIME TOTAL    / h     =   " << total << XMPI::endl << XMPI::endl; 
            XMPI::cout << ss.str();
        }
        // learn wisdom
        mDomain->learnWisdom(tstep - 1);
        
        // assemble phase 2: wait + extract 
        mDomain->assembleStiff(1);
    }
    ////////////////////////// loop //////////////////////////
    mDomain->dumpLeft();
    mDomain->dumpWisdom();
    double elapsed = timer.elapsed() * sec2h;
    if (verbose) {
        XMPI::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" << XMPI::endl;
        XMPI::cout << "TTTTTTTTT  NEWMARK TIME LOOP FINISHES  TTTTTTTTT" << XMPI::endl;
        XMPI::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" << XMPI::endl << XMPI::endl;
        XMPI::cout << "SIMULATION TIME / sec   =   " << t << XMPI::endl;
        XMPI::cout << "TOTAL STEPS DONE        =   " << maxStep << " / " <<  maxStep << XMPI::endl;
        XMPI::cout << "WALLTIME ELAPSED / h    =   " << elapsed << XMPI::endl << XMPI::endl;
        XMPI::cout << mDomain->reportCost();
    }
}



