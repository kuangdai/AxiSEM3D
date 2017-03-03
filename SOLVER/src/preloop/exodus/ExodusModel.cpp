// ExodusModel.cpp
// created by Kuangdai on 2-May-2016 
// read Exodus mesh file

// Adapted from salvus project 
// ExodusModel.cpp by Michael Afanasiev


#include "ExodusModel.h"
#include <sstream>
#include <iomanip>
#include <cmath>

#include "XMPI.h"
#include "Mapping.h"
#include <boost/algorithm/string.hpp>

#include "AttBuilder.h"
#include "XMath.h"

#include "XTimer.h"

extern "C" {
    #include "hdf5.h"
};
#include "H5Reader.h"


ExodusModel::ExodusModel(const std::string &fileName): mExodusFileName(fileName) {
    std::vector<std::string> substrs;
    boost::split(substrs, mExodusFileName, boost::is_any_of("/"), boost::token_compress_on);
    mExodusTitle = substrs[substrs.size() - 1];
    boost::trim_if(mExodusTitle, boost::is_any_of("\t "));
}

void ExodusModel::initialize() {
    XTimer::begin("Read Exodus", 1);
    if (XMPI::root()) readRawData();
    XTimer::end("Read Exodus", 1);
    
    XTimer::begin("Bcast Exodus", 1);
    bcastRawData();
    XTimer::end("Bcast Exodus", 1);
    
    XTimer::begin("Process Exodus", 1);
    formStructured();
    finishReading();
    XTimer::end("Process Exodus", 1);
}

void ExodusModel::readRawData() {
    /////////// open file
    hid_t fid = H5Fopen(mExodusFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fid < 0) throw std::runtime_error("ExodusModel::readRawData || "
        "Error opening exodus model file: ||" + mExodusFileName);
    
    /////////// file attr
    H5Reader::getAttribute(fid, "version", &mExodusVersion);
    
    /////////// numbers
    mNumNodes = H5Reader::getDimension1D(fid, "num_nodes");
    mNumQuads = H5Reader::getDimension1D(fid, "num_elem");
    
    /////////// global var
    int row, col;
    H5Reader::getStringData(fid, "name_glo_var", mRawNumGlbVar, mRawLenGlbVarName, mRawGlbVarName);
    H5Reader::getDoubleData(fid, "vals_glo_var", row, col, mRawGlbVarVals);
    
    /////////// global rec
    H5Reader::getStringData(fid, "info_records", mRawNumGlbRec, mRawLenGlbRec, mRawGlbRec);
    
    /////////// connectivity
    H5Reader::getIntData(fid, "connect1", row, col, mRawConnect);
    
    /////////// coords
    H5Reader::getDoubleData(fid, "coordx", row, col, mRawNodalS);
    H5Reader::getDoubleData(fid, "coordy", row, col, mRawNodalZ);

    /////////// elemental
    H5Reader::getStringData(fid, "name_elem_var", mRawNumEleVar, mRawLenEleVarName, mRawEleVarName);
    mRawEleVarVals = new double * [mRawNumEleVar];
    mRawEleVarVals[0] = new double[mRawNumEleVar * mNumQuads];
    for (int i = 1; i < mRawNumEleVar; i++) mRawEleVarVals[i] = mRawEleVarVals[0] + i * mNumQuads;
    for (int i = 0; i < mRawNumEleVar; i++) {
        std::stringstream ss;
        ss << "vals_elem_var" << i + 1 << "eb1";
        H5Reader::getDoubleData(fid, ss.str().c_str(), row, col, mRawEleVarVals[i], false);
    }
    
    /////////// side sets
    H5Reader::getStringData(fid, "ss_names", mRawNumSS, mRawLenSSName, mRawNameSS);
    mRawNumPairsSS = new int[mRawNumSS];
    mRawMaxPair = -1;
    for (int i = 0; i < mRawNumSS; i++) {
        std::stringstream ss;
        ss << "num_side_ss" << i + 1;
        mRawNumPairsSS[i] = H5Reader::getDimension1D(fid, ss.str().c_str());
        mRawMaxPair = std::max(mRawMaxPair, mRawNumPairsSS[i]);
    }
    mRawElemSS = new int * [mRawNumSS];
    mRawSideSS = new int * [mRawNumSS];
    mRawElemSS[0] = new int[mRawNumSS * mRawMaxPair];
    mRawSideSS[0] = new int[mRawNumSS * mRawMaxPair];
    for (int i = 0; i < mRawNumSS * mRawMaxPair; i++) {
        mRawElemSS[0][i] = -1;
        mRawSideSS[0][i] = -1;
    }
    for (int i = 1; i < mRawNumSS; i++) {
        mRawElemSS[i] = mRawElemSS[0] + i * mRawMaxPair;
        mRawSideSS[i] = mRawSideSS[0] + i * mRawMaxPair;
    }
    for (int i = 0; i < mRawNumSS; i++) {
        std::stringstream sse, sss;
        sse << "elem_ss" << i + 1;
        sss << "side_ss" << i + 1;
        H5Reader::getIntData(fid, sse.str().c_str(), row, col, mRawElemSS[i], false);
        H5Reader::getIntData(fid, sss.str().c_str(), row, col, mRawSideSS[i], false);
    }
    
    /////////// ellipticity
    H5Reader::getDoubleData(fid, "ellipticity", row, mRawEllipCol, mRawEllipData);
    
    /////////// close file
    H5Reader::hdf5Error(H5Fclose(fid), "H5Dclose");
}

void ExodusModel::bcastRawData() {
    /////////// file attr
    XMPI::bcast(mExodusVersion);
    
    /////////// numbers
    XMPI::bcast(mNumNodes);
    XMPI::bcast(mNumQuads);
    
    /////////// global var
    XMPI::bcast(mRawNumGlbVar);
    XMPI::bcast(mRawLenGlbVarName);
    
    XMPI::bcast_alloc(mRawGlbVarName, mRawNumGlbVar * mRawLenGlbVarName);
    XMPI::bcast_alloc(mRawGlbVarVals, mRawNumGlbVar);
    
    /////////// global rec
    XMPI::bcast(mRawNumGlbRec);
    XMPI::bcast(mRawLenGlbRec);
    XMPI::bcast_alloc(mRawGlbRec, mRawNumGlbRec * mRawLenGlbRec);
    
    /////////// connectivity
    XMPI::bcast_alloc(mRawConnect, mNumQuads * 4);
    
    /////////// coords
    XMPI::bcast_alloc(mRawNodalS, mNumNodes);
    XMPI::bcast_alloc(mRawNodalZ, mNumNodes);
    
    /////////// elemental
    XMPI::bcast(mRawNumEleVar);
    XMPI::bcast(mRawLenEleVarName);
    XMPI::bcast_alloc(mRawEleVarName, mRawNumEleVar * mRawLenEleVarName);
    if (!XMPI::root()) {
        mRawEleVarVals = new double * [mRawNumEleVar];
        mRawEleVarVals[0] = new double[mRawNumEleVar * mNumQuads];
        for (int i = 1; i < mRawNumEleVar; i++) mRawEleVarVals[i] = mRawEleVarVals[0] + i * mNumQuads;
    }
    XMPI::bcast(mRawEleVarVals[0], mRawNumEleVar * mNumQuads);
    
    /////////// side sets
    XMPI::bcast(mRawNumSS);
    XMPI::bcast(mRawLenSSName);
    XMPI::bcast(mRawMaxPair);
    XMPI::bcast_alloc(mRawNameSS, mRawNumSS * mRawLenSSName);
    XMPI::bcast_alloc(mRawNumPairsSS, mRawNumSS);
    if (!XMPI::root()) {
        mRawElemSS = new int * [mRawNumSS];
        mRawSideSS = new int * [mRawNumSS];
        mRawElemSS[0] = new int[mRawNumSS * mRawMaxPair];
        mRawSideSS[0] = new int[mRawNumSS * mRawMaxPair];
        for (int i = 1; i < mRawNumSS; i++) {
            mRawElemSS[i] = mRawElemSS[0] + i * mRawMaxPair;
            mRawSideSS[i] = mRawSideSS[0] + i * mRawMaxPair;
        }
    }
    XMPI::bcast(mRawElemSS[0], mRawNumSS * mRawMaxPair);
    XMPI::bcast(mRawSideSS[0], mRawNumSS * mRawMaxPair);
    
    /////////// ellipticity
    XMPI::bcast(mRawEllipCol);
    XMPI::bcast_alloc(mRawEllipData, mRawEllipCol * 2);
}

void ExodusModel::formStructured() {
    /////////// global var
    std::string all(mRawGlbVarName);
    for (int i = 0; i < mRawNumGlbVar; i++) {
        std::string varName = all.substr(i * mRawLenGlbVarName, mRawLenGlbVarName);
        boost::trim_if(varName, boost::is_any_of("\t "));
        if (varName == "dt") varName = "dt (nPol = 1)";
        mGlobalVariables.insert(std::pair<std::string, double>(varName, mRawGlbVarVals[i]));
    }
    delete [] mRawGlbVarName;
    delete [] mRawGlbVarVals;
    
    /////////// global rec
    all = std::string(mRawGlbRec);
    std::vector<std::string> included = {"crdsys", "model"};
    for (int i = 0; i < mRawNumGlbRec; i++) {
        std::string recstr = all.substr(i * mRawLenGlbRec, mRawLenGlbRec);
        std::vector<std::string> substrs;
        boost::trim_if(recstr, boost::is_any_of("\t "));
        boost::split(substrs, recstr, boost::is_any_of("="), boost::token_compress_on);
        if (std::find(included.begin(), included.end(), boost::trim_copy(substrs[0])) != included.end()) {
            mGlobalRecords.insert(std::pair<std::string, std::string>(
                boost::trim_copy(substrs[0]), boost::trim_copy(substrs[1])));
            if (substrs[0].find("crdsys") != std::string::npos) 
                mCartesian = (substrs[1].find("cartesian") != std::string::npos);
        }
    }
    delete [] mRawGlbRec;
    // name of axis and surface sets
    mSSNameAxis = mCartesian ? "x0" : "t0";
    mSSNameSurface = mCartesian ? "y1" : "r1";
    
    /////////// connectivity
    mConnectivity.resize(mNumQuads);
    for (int i = 0; i < mNumQuads; i++) {
        mConnectivity[i][0] = mRawConnect[i * 4 + 0] - 1;
        mConnectivity[i][1] = mRawConnect[i * 4 + 1] - 1;
        mConnectivity[i][2] = mRawConnect[i * 4 + 2] - 1;
        mConnectivity[i][3] = mRawConnect[i * 4 + 3] - 1;
    }
    delete [] mRawConnect;
    
    /////////// coords
    mNodalS.resize(mNumNodes);
    mNodalZ.resize(mNumNodes);
    for (int i = 0; i < mNumNodes; i++) {
        mNodalS[i] = mRawNodalS[i];
        mNodalZ[i] = mRawNodalZ[i];
    }
    // NOTE: we temporarily treat Cartesian meshes as special cases of spherical meshes
    //       by means of moving it to the "north pole". The introduced global 
    //       curvature should be ignorable, or the problem itself is ill-defined
    //       as a local problem.
    if (mCartesian) {
        double R_EARTH = mGlobalVariables.at("radius");
        double maxz = -1.;
        for (int i = 0; i < mNumNodes; i++) 
            maxz = std::max(maxz, mNodalZ[i]);
        for (int i = 0; i < mNumNodes; i++) 
            mNodalZ[i] += R_EARTH - maxz;
    }
    delete [] mRawNodalS;
    delete [] mRawNodalZ;
    
    /////////// elemental
    all = std::string(mRawEleVarName);
    for (int i = 0; i < mRawNumEleVar; i++) {
        std::string varName = all.substr(i * mRawLenEleVarName, mRawLenEleVarName);
        boost::trim_if(varName, boost::is_any_of("\t "));
        std::vector<double> buffer(mNumQuads);
        for (int j = 0; j < mNumQuads; j++) buffer[j] = mRawEleVarVals[i][j];
        mElementalVariables.insert(std::pair<std::string, std::vector<double>>
                    (std::string(varName), buffer));
    }
    delete [] mRawEleVarName;
    delete [] mRawEleVarVals[0];
    delete [] mRawEleVarVals;
    
    /////////// side sets 
    all = std::string(mRawNameSS);
    for (int i = 0; i < mRawNumSS; i++) {
        std::string varName = all.substr(i * mRawLenSSName, mRawLenSSName);
        boost::trim_if(varName, boost::is_any_of("\t "));
        std::vector<int> ss(mNumQuads, -1);
        for (int j = 0; j < mRawNumPairsSS[i]; j++) ss[mRawElemSS[i][j] - 1] = mRawSideSS[i][j] - 1;
        mSideSets.insert(std::pair<std::string, std::vector<int>>(std::string(varName), ss));
    }
    delete [] mRawNameSS;
    delete [] mRawNumPairsSS;
    delete [] mRawElemSS[0];
    delete [] mRawSideSS[0];
    delete [] mRawElemSS;
    delete [] mRawSideSS;
    
    /////////// ellipticity
    for (int i = 0; i < mRawEllipCol; i++) {
        mEllipKnots.push_back(mRawEllipData[i]);
        mEllipCoeffs.push_back(mRawEllipData[i + mRawEllipCol]);
    }
    delete [] mRawEllipData;
}

void ExodusModel::finishReading() {
    XTimer::begin("Process Exodus DistTol", 2);
    // distance tolerance
    double distTol = 1e100;
    for (int i = 0; i < mNumQuads; i++) {
        if (i % XMPI::nproc() != XMPI::rank()) continue;
        double s0 = mNodalS[mConnectivity[i][0]];
        double z0 = mNodalZ[mConnectivity[i][0]];
        double s1 = mNodalS[mConnectivity[i][1]];
        double z1 = mNodalZ[mConnectivity[i][1]];
        double s2 = mNodalS[mConnectivity[i][2]];
        double z2 = mNodalZ[mConnectivity[i][2]];
        double s3 = mNodalS[mConnectivity[i][3]];
        double z3 = mNodalZ[mConnectivity[i][3]];
        double dist0 = sqrt((s0 - s1) * (s0 - s1) + (z0 - z1) * (z0 - z1)) / 1000.;
        double dist1 = sqrt((s1 - s2) * (s1 - s2) + (z1 - z2) * (z1 - z2)) / 1000.;
        double dist2 = sqrt((s2 - s3) * (s2 - s3) + (z2 - z3) * (z2 - z3)) / 1000.;
        double dist3 = sqrt((s3 - s0) * (s3 - s0) + (z3 - z0) * (z3 - z0)) / 1000.;
        distTol = std::min({dist0, dist1, dist2, dist3, distTol});
    }
    mDistTolerance = XMPI::min(distTol);
    XTimer::end("Process Exodus DistTol", 2);
    
    // surface radius
    XTimer::begin("Process Exodus ROuter", 2);
    double router = -1.;
    for (int i = 0; i < mNumNodes; i++) {
        if (i % XMPI::nproc() != XMPI::rank()) continue;
        router = std::max(router, mNodalZ[i]);
    } 
    mROuter = XMPI::max(router);
    XTimer::end("Process Exodus ROuter", 2);
        
    // average gll spacing
    XTimer::begin("Process Exodus GLL-Spacing", 2);
    std::vector<std::vector<int>> refElem(mNumNodes, std::vector<int>());
    for (int i = 0; i < mNumQuads; i++) {
        refElem[mConnectivity[i][0]].push_back(i);
        refElem[mConnectivity[i][1]].push_back(i);
        refElem[mConnectivity[i][2]].push_back(i);
        refElem[mConnectivity[i][3]].push_back(i);
    }
    mAveGLLSpacing = std::vector<double>(mNumNodes, 0.);
    for (int i = 0; i < mNumNodes; i++) {
        if (i % XMPI::nproc() != XMPI::rank()) continue;
        for (int j = 0; j < refElem[i].size(); j++) {
            int ielem = refElem[i][j];
            double s0 = mNodalS[mConnectivity[ielem][0]];
            double z0 = mNodalZ[mConnectivity[ielem][0]];
            double s1 = mNodalS[mConnectivity[ielem][1]];
            double z1 = mNodalZ[mConnectivity[ielem][1]];
            double s2 = mNodalS[mConnectivity[ielem][2]];
            double z2 = mNodalZ[mConnectivity[ielem][2]];
            double s3 = mNodalS[mConnectivity[ielem][3]];
            double z3 = mNodalZ[mConnectivity[ielem][3]];
            double dist0 = sqrt((s0 - s1) * (s0 - s1) + (z0 - z1) * (z0 - z1));
            double dist1 = sqrt((s1 - s2) * (s1 - s2) + (z1 - z2) * (z1 - z2));
            double dist2 = sqrt((s2 - s3) * (s2 - s3) + (z2 - z3) * (z2 - z3));
            double dist3 = sqrt((s3 - s0) * (s3 - s0) + (z3 - z0) * (z3 - z0));
            mAveGLLSpacing[i] += (dist0 + dist1 + dist2 + dist3) / 4. / nPol / refElem[i].size();
        }
    } 
    XMPI::sumVector(mAveGLLSpacing);
    XTimer::end("Process Exodus GLL-Spacing", 2);
    
    // rotate nodes of axial elements such that side 3 is on axis
    XTimer::begin("Process Exodus Axis", 2);
    for (int axialQuad = 0; axialQuad < mNumQuads; axialQuad++) {
        // loop over t0
        int axialSide = getSideAxis(axialQuad);
        if (axialSide == 3 || axialSide == -1) continue;

        // connectivity
        std::array<int, 4> con = mConnectivity[axialQuad];
        for (int j = 0; j < 4; j++) 
            mConnectivity[axialQuad][j] = con[Mapping::period0123(j + axialSide - 3)];
        
        
        // elemental fields
        for (auto it = mElementalVariables.begin(); it != mElementalVariables.end(); it++) {
            std::string vname = it->first;
            if (vname.substr(vname.length() - 2, 2) == std::string("_0")) {
                vname = vname.substr(0, vname.length() - 2);
                std::array<double, 4> v_old;
                v_old[0] = mElementalVariables.at(vname + "_0")[axialQuad];
                v_old[1] = mElementalVariables.at(vname + "_1")[axialQuad];
                v_old[2] = mElementalVariables.at(vname + "_2")[axialQuad];
                v_old[3] = mElementalVariables.at(vname + "_3")[axialQuad];
                mElementalVariables.at(vname + "_0")[axialQuad] = v_old[Mapping::period0123(0 + axialSide - 3)];
                mElementalVariables.at(vname + "_1")[axialQuad] = v_old[Mapping::period0123(1 + axialSide - 3)];
                mElementalVariables.at(vname + "_2")[axialQuad] = v_old[Mapping::period0123(2 + axialSide - 3)];
                mElementalVariables.at(vname + "_3")[axialQuad] = v_old[Mapping::period0123(3 + axialSide - 3)];
            }
        }
        
        // side sets
        for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
            if (it->second[axialQuad] != -1) {
                it->second[axialQuad] = Mapping::period0123(it->second[axialQuad] - axialSide + 3);
            } 
        } 
        
        // done
        // mSideSets.at(mSSNameAxis)[axialQuad] = 3;
    }
    XTimer::end("Process Exodus Axis", 2);
    
    // find elements that are not axial but neighboring axial elements
    XTimer::begin("Process Exodus Vicinal", 2);
    std::array<int, 4> data = {-1, -1, -1, -1};
    mVicinalAxis = std::vector<std::array<int, 4>>(mNumQuads, data);
    // first find near-axis nodes and axial quads
    std::vector<bool> nodeNearAxis(mNumNodes, false);
    std::vector<bool> quadOnAxis(mNumQuads, false);
    for (int axialQuad = 0; axialQuad < mNumQuads; axialQuad++) {
        int axialSide = getSideAxis(axialQuad);
        if (axialSide == -1) continue;
        nodeNearAxis[mConnectivity[axialQuad][0]] = true;
        nodeNearAxis[mConnectivity[axialQuad][1]] = true;
        nodeNearAxis[mConnectivity[axialQuad][2]] = true;
        nodeNearAxis[mConnectivity[axialQuad][3]] = true;
        quadOnAxis[axialQuad] = true;
    }
    // loop over quads
    for (int iquad = 0; iquad < mNumQuads; iquad++) {
        if (quadOnAxis[iquad]) continue;
        for (int j = 0; j < 4; j++) {
            int nTag = mConnectivity[iquad][j];
            if (nodeNearAxis[nTag]) mVicinalAxis[iquad][j] = j;
        }
    }
    XTimer::end("Process Exodus Vicinal", 2);

    // check if ocean presents in mesh
    XTimer::begin("Process Exodus Check Ocean", 2);
    std::string strVs = isIsotropic() ? "VS_0" : "VSV_0";
    for (int iQuad = 0; iQuad < mNumQuads; iQuad++) {
        if (iQuad % XMPI::nproc() != XMPI::rank()) continue;
        if (getSideSurface(iQuad) == -1) continue;
        double vs = mElementalVariables.at(strVs)[iQuad];
        if (vs < tinyDouble) throw std::runtime_error("ExodusModel::finishReading || "
            "Ocean is detected in mesh. By far, realistic ocean is not implemented in AxiSEM3D. ||"
            "Use non-ocean models in the Mesher and add ocean load in inparam.basic.");
    }
    XTimer::end("Process Exodus Check Ocean", 2);
}

std::string ExodusModel::verbose() const {
    std::stringstream ss;
    ss << "\n======================= Exodus Model =======================" << std::endl;
    ss << "  Overview__________________________________________________" << std::endl;
    ss << "    Exodus Title      =   " << mExodusTitle << std::endl;
    ss << "    Mesh CS Type      =   " << (mCartesian ? "Cartesian" : "Spherical") << std::endl;
    ss << "    Number of Nodes   =   " << mNumNodes << std::endl;
    ss << "    Number of Quads   =   " << mNumQuads << std::endl;
    ss << "  Global Variables__________________________________________" << std::endl;
    int widthname = -1;
    for (auto it = mGlobalVariables.begin(); it != mGlobalVariables.end(); it++) 
        widthname = std::max(widthname, (int)(it->first.length()));
    for (auto it = mGlobalRecords.begin(); it != mGlobalRecords.end(); it++) 
        widthname = std::max(widthname, (int)(it->first.length()));    
    for (auto it = mGlobalVariables.begin(); it != mGlobalVariables.end(); it++) 
        ss << "    " << std::setw(widthname) << it->first << "   =   " << it->second << std::endl;
    for (auto it = mGlobalRecords.begin(); it != mGlobalRecords.end(); it++) 
        ss << "    " << std::setw(widthname) << it->first << "   =   " << it->second << std::endl;
    ss << "  Connectivity______________________________________________" << std::endl;
    int width = (int)std::log10(std::max(mNumQuads, mNumNodes)) + 1;
    ss << "    " << std::setw(width) << 0 << ": ";
    for (int j = 0; j < 4; j++) ss << std::setw(width) << mConnectivity[0][j] << " ";
    ss << std::endl << "    " << std::setw(width) << "..." << std::endl;
    ss << "    " << std::setw(width) << mNumQuads - 1 << ": ";
    for (int j = 0; j < 4; j++) ss << std::setw(width) << mConnectivity[mNumQuads - 1][j] << " ";
    ss << std::endl;
    ss << "  Coordinates_______________________________________________" << std::endl;
    ss << "    " << std::setw(width) << 0 << ": ";
    ss << std::setw(14) << mNodalS[0] << std::setw(14) << mNodalZ[0] << std::endl;
    ss<< "    " << std::setw(width) << "..." << std::endl;
    ss << "    " << std::setw(width) << mNumNodes - 1 << ": ";
    ss << std::setw(14) << mNodalS[mNumNodes - 1] << std::setw(14) << mNodalZ[mNumNodes - 1] << std::endl;
    ss << "  Elemental Variables_______________________________________" << std::endl;
    widthname = -1;
    for (auto it = mElementalVariables.begin(); it != mElementalVariables.end(); it++) 
        widthname = std::max(widthname, (int)(it->first.length()));
    for (auto it = mElementalVariables.begin(); it != mElementalVariables.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << ": ";
        ss << std::setw(14) << it->second[0] << ", ..., ";
        ss << std::setw(14) << it->second[mNumQuads - 1] << std::endl;        
    }
    ss << "  Side Sets_________________________________________________" << std::endl;
    widthname = -1;
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) 
        widthname = std::max(widthname, (int)(it->first.length()));
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << ":   ";
        int pair = 0;
        for (int q = 0; q < mNumQuads; q++) pair += (int)(it->second[q] >= 0);
        ss << pair << " paris" << std::endl;
        // ss << "  (" << std::setw(width) << it->second[0][0] << ",";
        // ss << std::setw(width) << it->second[0][1] << ")" << ", ..., ";
        // ss << "(" << std::setw(width) << it->second[it->second.size() - 1][0] << ",";
        // ss << std::setw(width) << it->second[it->second.size() - 1][1] << ")" << std::endl;        
    }
    ss << "  Miscellaneous_____________________________________________" << std::endl;
    ss << "    Distance Tolerance / m   =   " << mDistTolerance << std::endl;
    ss << "    Surface Radius / m       =   " << mROuter << std::endl;
    if (!mCartesian) {
        ss << "  External__________________________________________________" << std::endl;
        ss << "    Num. Ellipticity Spline Knots   =   " << mEllipKnots.size() << std::endl;
    }
    ss << "======================= Exodus Model =======================\n" << std::endl;
    return ss.str();
}

#include "Parameters.h"
#include "global.h"
void ExodusModel::buildInparam(ExodusModel *&exModel, const Parameters &par, 
    AttParameters *&attPar, int verbose) {
    if (exModel) delete exModel;
    std::string exfile = par.getValue<std::string>("MODEL_1D_EXODUS_MESH_FILE");
    exfile = Parameters::sInputDirectory + "/" + exfile;
    exModel = new ExodusModel(exfile);
    exModel->initialize();
    if (verbose) XMPI::cout << exModel->verbose();
    
    // form attenuation parameters
    int nr_lin_solids = (int)exModel->mGlobalVariables.at("nr_lin_solids");
    double f_min = exModel->mGlobalVariables.at("f_min");
    double f_max = exModel->mGlobalVariables.at("f_max");
    double f_ref = exModel->mGlobalVariables.at("f_ref");
    RDColX w(nr_lin_solids), y(nr_lin_solids);
    for (int i = 0; i < nr_lin_solids; i++) {
        std::stringstream sw, sy;
        sw << "w_" << i;
        sy << "y_" << i;
        w(i) = exModel->mGlobalVariables.at(sw.str());
        y(i) = exModel->mGlobalVariables.at(sy.str());
    }
    if (attPar) delete attPar;
    attPar = new AttParameters(nr_lin_solids, f_min, f_max, f_ref, w, y);
    
    // ellipticity
    if (exModel->isCartesian()) return;
    std::string emode = par.getValue<std::string>("MODEL_3D_ELLIPTICITY_MODE");
    if (boost::iequals(emode, "off")) {
        // no ellipticity
        XMath::setEllipticity(0., exModel->mROuter, std::vector<double>(), std::vector<double>());
    } else {
        double inv_f = par.getValue<double>("MODEL_3D_ELLIPTICITY_INVF");
        if (inv_f <= 0.) throw std::runtime_error("ExodusModel::buildInparam || Invalid flattening."); 
        XMath::setEllipticity(1. / inv_f, exModel->mROuter, exModel->mEllipKnots, exModel->mEllipCoeffs);
    }
}

