// axisem.h
// created by Kuangdai on 26-Mar-2016 
// main of AxiSEM3D

#pragma once

#include <vector>
#include "XMPI.h"

// preloop 
#include "SpectralConstants.h"
#include "Parameters.h"
#include "ExodusModel.h"
#include "NrField.h"
#include "Volumetric3D.h"
#include "Geometric3D.h"
#include "OceanLoad3D.h"
#include "Source.h"
#include "Mesh.h"
#include "AttBuilder.h"
#include "STF.h"
#include "ReceiverCollection.h"

// solver
#include "Domain.h"
#include "Newmark.h"

struct PreloopVariables {
    Parameters *mParameters = 0;
    AttParameters *mAttParameters = 0;
    ExodusModel *mExodusModel = 0;
    NrField *mNrField = 0;
    std::vector<Volumetric3D *> mVolumetric3D;
    std::vector<Geometric3D *> mGeometric3D;
    OceanLoad3D *mOceanLoad3D = 0;
    Source *mSource = 0;
    Mesh *mMesh = 0;
    AttBuilder *mAttBuilder = 0;
    STF *mSTF = 0;
    ReceiverCollection *mReceivers = 0;
    
    // finalizer 
    void finalize() {
        if (mParameters) {delete mParameters; mParameters = 0;}
        if (mAttParameters) {delete mAttParameters; mAttParameters = 0;}
        if (mExodusModel) {delete mExodusModel; mExodusModel = 0;}
        if (mNrField) {delete mNrField; mNrField = 0;}
        for (const auto &m: mVolumetric3D) delete m;
        mVolumetric3D.clear();
        for (const auto &m: mGeometric3D) delete m;
        mGeometric3D.clear();
        if (mOceanLoad3D) {delete mOceanLoad3D; mOceanLoad3D = 0;}
        if (mSource) {delete mSource; mSource = 0;}
        if (mMesh) {delete mMesh; mMesh = 0;}
        if (mAttBuilder) {delete mAttBuilder; mAttBuilder = 0;}
        if (mSTF) {delete mSTF; mSTF = 0;}
        if (mReceivers) {delete mReceivers; mReceivers = 0;}
    };
};

struct SolverVariables {
    Domain *mDomain = 0;
    Newmark *mNewmark = 0;
    
    // finalizer
    void finalize() {
        if (mDomain) {delete mDomain; mDomain = 0;}
        if (mNewmark) {delete mNewmark; mNewmark = 0;}
    };
};

//////////////////////////////// functons ////////////////////////////////
int axisem_main(int argc, char *argv[]);
void initializeSolverStatic(int maxNr, bool disableWisdomFFTW);
void finalizeSolverStatic();



