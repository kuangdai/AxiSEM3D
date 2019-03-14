// ExodusModel.cpp
// created by Kuangdai on 2-May-2016 
// read Exodus mesh file

// Adapted from salvus project 
// ExodusModel.cpp by Michael Afanasiev


#include "ExodusModel.h"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat>

#include "XMPI.h"
#include "Mapping.h"
#include <boost/algorithm/string.hpp>

#include "AttBuilder.h"
#include "Geodesy.h"
#include "Parameters.h"
#include "MultilevelTimer.h"

#include "NetCDF_Reader.h"

ExodusModel::ExodusModel(const std::string &fileName): mExodusFileName(fileName) {
    std::vector<std::string> substrs = Parameters::splitString(mExodusFileName, "/");
    mExodusTitle = substrs[substrs.size() - 1];
    boost::trim_if(mExodusTitle, boost::is_any_of("\t "));
}

void ExodusModel::initialize() {
    MultilevelTimer::begin("Read Exodus", 1);
    if (XMPI::root()) {
        readRawData();
    }
    MultilevelTimer::end("Read Exodus", 1);
    
    MultilevelTimer::begin("Bcast Exodus", 1);
    bcastRawData();
    MultilevelTimer::end("Bcast Exodus", 1);
    
    MultilevelTimer::begin("Process Exodus", 1);
    formStructured();
    formAuxiliary();
    MultilevelTimer::end("Process Exodus", 1);
}

void ExodusModel::readRawData() {
    // open file
    NetCDF_Reader reader;
    reader.open(mExodusFileName);
    RDMatXX dbuffer;
    
    // global
    reader.readString("name_glo_var", mGlobalVariableNames);
    reader.read1D("vals_glo_var", mGlobalVariableValues);
    reader.readString("info_records", mGlobalRecordsRaw);
    
    // connectivity and coords
    reader.read2D("connect1", mConnectivity);
    mConnectivity.array() -= 1;
    reader.read1D("coordx", mNodalS);
    reader.read1D("coordy", mNodalZ);
        
    // elemental variables
    reader.readString("name_elem_var", mElementalVariableNames);
    mElementalVariableValues = RDMatXX::Zero(getNumQuads(), mElementalVariableNames.size());
    for (int i = 0; i < mElementalVariableNames.size(); i++) {
        std::stringstream ss;
        ss << "vals_elem_var" << i + 1 << "eb1";
        reader.read2D(ss.str(), dbuffer);
        mElementalVariableValues.col(i) = dbuffer.transpose();
    }    
    
    // side sets
    reader.readString("ss_names", mSideSetNames);
    mSideSetValues = IMatXX::Zero(getNumQuads(), mSideSetNames.size());
    for (int i = 0; i < mSideSetNames.size(); i++) {
        std::stringstream sse, sss;
        IColX elems, sides;
        // elem
        sse << "elem_ss" << i + 1;
        reader.read1D(sse.str(), elems);
        // side
        sss << "side_ss" << i + 1;
        reader.read1D(sss.str(), sides);
        // processing
        IColX values = IColX::Constant(getNumQuads(), -1);
        for (int j = 0; j < sides.rows(); j++) {
            values(elems(j) - 1) = sides(j) - 1;
        }
        mSideSetValues.col(i) = values;
    }
    
    // do not read ellipticity for cartesian
    bool cartesian = false;
    for (int i = 0; i < mGlobalRecordsRaw.size(); i++) {
        std::vector<std::string> substrs = Parameters::splitString(mGlobalRecordsRaw[i], "=");
        if (boost::iequals(boost::trim_copy(substrs[0]), "crdsys") && 
            boost::iequals(boost::trim_copy(substrs[1]), "cartesian")) {
            cartesian = true;
            break;
        }
    }
    if (cartesian) {
        mEllipKnots = RDColX::Zero(0);
        mEllipCoeffs = RDColX::Zero(0);
        // set radius to PREM 
        mGlobalVariables.insert(std::pair<std::string, double>("radius", 6371e3));
    } else {
        // ellipticity
        reader.read2D("ellipticity", dbuffer);
        mEllipKnots = dbuffer.row(0).transpose();
        mEllipCoeffs = dbuffer.row(1).transpose();
    }
    
    // close file
    reader.close();
}

void ExodusModel::bcastRawData() {
    XMPI::bcast(mGlobalVariableNames);
    XMPI::bcastEigen(mGlobalVariableValues);
    XMPI::bcast(mGlobalRecordsRaw);
    
    XMPI::bcastEigen(mConnectivity);
    XMPI::bcastEigen(mNodalS);
    XMPI::bcastEigen(mNodalZ);
    
    XMPI::bcast(mElementalVariableNames);
    XMPI::bcastEigen(mElementalVariableValues);
    
    XMPI::bcast(mSideSetNames);
    XMPI::bcastEigen(mSideSetValues);
    
    XMPI::bcastEigen(mEllipKnots);
    XMPI::bcastEigen(mEllipCoeffs);
}

void ExodusModel::formStructured() {
    // global variables
    for (int i = 0; i < mGlobalVariableNames.size(); i++) {
        std::string varName = mGlobalVariableNames[i];
        if (varName == "dt") {
            varName = "dt (nPol = 1)";
        }
        mGlobalVariables.insert(std::pair<std::string, double>(varName, mGlobalVariableValues(i)));
    }
    
    // global records
    std::vector<std::string> included = {"crdsys", "model"};
    for (int i = 0; i < mGlobalRecordsRaw.size(); i++) {
        std::vector<std::string> substrs = Parameters::splitString(mGlobalRecordsRaw[i], "=");
        if (std::find(included.begin(), included.end(), boost::trim_copy(substrs[0])) != included.end()) {
            mGlobalRecords.insert(std::pair<std::string, std::string>(
                boost::trim_copy(substrs[0]), boost::trim_copy(substrs[1])));
        }
    }
    
    // elemental variables
    for (int i = 0; i < mElementalVariableNames.size(); i++) {
        mElementalVariables.insert(std::pair<std::string, RDColX>(mElementalVariableNames[i], 
            mElementalVariableValues.col(i)));
    }
    
    // side sets
    for (int i = 0; i < mSideSetNames.size(); i++) {
        mSideSets.insert(std::pair<std::string, IColX>(mSideSetNames[i], 
            mSideSetValues.col(i)));
    }
    
    // name of axis and surface sets
    std::string sphere_axis_name = "t0";
    if (!isCartesian()) {
        if (mSideSets.find("t1") != mSideSets.end()) {
            sphere_axis_name = "t1";
        }
    }
    
    mSSNameAxis = isCartesian() ? "x0" : sphere_axis_name;
    mSSNameSurface = isCartesian() ? "y1" : "r1";
    
    // NOTE: we temporarily treat Cartesian meshes as special cases of spherical meshes
    //       by means of moving it to the "north pole". The introduced global 
    //       curvature should be ignorable, or the problem itself is ill-defined
    //       as a local problem.
    if (isCartesian()) {
        double R_EARTH = getROuter();
        double maxz = mNodalZ.maxCoeff();
        mNodalZ.array() += R_EARTH - maxz; 
    }
}

void ExodusModel::formAuxiliary() {
    MultilevelTimer::begin("Process Exodus DistTol", 2);
    // distance tolerance
    double distTol = DBL_MAX;
    for (int i = 0; i < getNumQuads(); i++) {
        if (i % XMPI::nproc() != XMPI::rank()) {
            continue;
        }
        double s0 = mNodalS(mConnectivity(i, 0));
        double z0 = mNodalZ(mConnectivity(i, 0));
        double s1 = mNodalS(mConnectivity(i, 1));
        double z1 = mNodalZ(mConnectivity(i, 1));
        double s2 = mNodalS(mConnectivity(i, 2));
        double z2 = mNodalZ(mConnectivity(i, 2));
        double s3 = mNodalS(mConnectivity(i, 3));
        double z3 = mNodalZ(mConnectivity(i, 3));
        double dist0 = sqrt((s0 - s1) * (s0 - s1) + (z0 - z1) * (z0 - z1)) / 1000.;
        double dist1 = sqrt((s1 - s2) * (s1 - s2) + (z1 - z2) * (z1 - z2)) / 1000.;
        double dist2 = sqrt((s2 - s3) * (s2 - s3) + (z2 - z3) * (z2 - z3)) / 1000.;
        double dist3 = sqrt((s3 - s0) * (s3 - s0) + (z3 - z0) * (z3 - z0)) / 1000.;
        distTol = std::min({dist0, dist1, dist2, dist3, distTol});
    }
    mDistTolerance = XMPI::min(distTol);
    MultilevelTimer::end("Process Exodus DistTol", 2);
    
    // average gll spacing
    MultilevelTimer::begin("Process Exodus GLL-Spacing", 2);
    std::vector<std::vector<int>> refElem(getNumNodes(), std::vector<int>());
    for (int i = 0; i < getNumQuads(); i++) {
        refElem[mConnectivity(i, 0)].push_back(i);
        refElem[mConnectivity(i, 1)].push_back(i);
        refElem[mConnectivity(i, 2)].push_back(i);
        refElem[mConnectivity(i, 3)].push_back(i);
    }
    mAveGLLSpacing = RDColX::Zero(getNumNodes());
    for (int i = 0; i < getNumNodes(); i++) {
        if (i % XMPI::nproc() != XMPI::rank()) {
            continue;
        }
        for (int j = 0; j < refElem[i].size(); j++) {
            int ielem = refElem[i][j];
            double s0 = mNodalS(mConnectivity(ielem, 0));
            double z0 = mNodalZ(mConnectivity(ielem, 0));
            double s1 = mNodalS(mConnectivity(ielem, 1));
            double z1 = mNodalZ(mConnectivity(ielem, 1));
            double s2 = mNodalS(mConnectivity(ielem, 2));
            double z2 = mNodalZ(mConnectivity(ielem, 2));
            double s3 = mNodalS(mConnectivity(ielem, 3));
            double z3 = mNodalZ(mConnectivity(ielem, 3));
            double dist0 = sqrt((s0 - s1) * (s0 - s1) + (z0 - z1) * (z0 - z1));
            double dist1 = sqrt((s1 - s2) * (s1 - s2) + (z1 - z2) * (z1 - z2));
            double dist2 = sqrt((s2 - s3) * (s2 - s3) + (z2 - z3) * (z2 - z3));
            double dist3 = sqrt((s3 - s0) * (s3 - s0) + (z3 - z0) * (z3 - z0));
            mAveGLLSpacing(i) += (dist0 + dist1 + dist2 + dist3) / 4. / nPol / refElem[i].size();
        }
    } 
    XMPI::sumEigenDouble(mAveGLLSpacing);
    MultilevelTimer::end("Process Exodus GLL-Spacing", 2);
    
    // rotate nodes of axial elements such that side 3 is on axis
    MultilevelTimer::begin("Process Exodus Axis", 2);
    for (int axialQuad = 0; axialQuad < getNumQuads(); axialQuad++) {
        // loop over t0
        int axialSide = getSideAxis(axialQuad);
        if (axialSide == 3 || axialSide == -1) {
            continue;
        }
    
        // connectivity
        IRow4 con = mConnectivity.row(axialQuad);
        for (int j = 0; j < 4; j++) {
            mConnectivity(axialQuad, j) = con(Mapping::period0123(j + axialSide - 3));
        }
        
        
        // elemental fields
        for (auto it = mElementalVariables.begin(); it != mElementalVariables.end(); it++) {
            std::string vname = it->first;
            if (vname.substr(vname.length() - 2, 2) == std::string("_0")) {
                vname = vname.substr(0, vname.length() - 2);
                RDRow4 v_old;
                v_old(0) = mElementalVariables.at(vname + "_0")(axialQuad);
                v_old(1) = mElementalVariables.at(vname + "_1")(axialQuad);
                v_old(2) = mElementalVariables.at(vname + "_2")(axialQuad);
                v_old(3) = mElementalVariables.at(vname + "_3")(axialQuad);
                mElementalVariables.at(vname + "_0")(axialQuad) = v_old(Mapping::period0123(0 + axialSide - 3));
                mElementalVariables.at(vname + "_1")(axialQuad) = v_old(Mapping::period0123(1 + axialSide - 3));
                mElementalVariables.at(vname + "_2")(axialQuad) = v_old(Mapping::period0123(2 + axialSide - 3));
                mElementalVariables.at(vname + "_3")(axialQuad) = v_old(Mapping::period0123(3 + axialSide - 3));
            }
        }
        
        // side sets
        for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
            if (it->second(axialQuad) != -1) {
                it->second(axialQuad) = Mapping::period0123(it->second(axialQuad) - axialSide + 3);
            } 
        } 
        
        // done
        // mSideSets.at(mSSNameAxis)(axialQuad) = 3;
    }
    MultilevelTimer::end("Process Exodus Axis", 2);
    
    // find elements that are not axial but neighboring axial elements
    MultilevelTimer::begin("Process Exodus Vicinal", 2);
    mVicinalAxis = IMatX4::Constant(getNumQuads(), 4, -1);
    // first find near-axis nodes and axial quads
    std::vector<bool> nodeNearAxis(getNumNodes(), false);
    std::vector<bool> quadOnAxis(getNumQuads(), false);
    for (int axialQuad = 0; axialQuad < getNumQuads(); axialQuad++) {
        int axialSide = getSideAxis(axialQuad);
        if (axialSide == -1) {
            continue;
        }
        nodeNearAxis[mConnectivity(axialQuad, 0)] = true;
        nodeNearAxis[mConnectivity(axialQuad, 1)] = true;
        nodeNearAxis[mConnectivity(axialQuad, 2)] = true;
        nodeNearAxis[mConnectivity(axialQuad, 3)] = true;
        quadOnAxis[axialQuad] = true;
    }
    // loop over quads
    for (int iquad = 0; iquad < getNumQuads(); iquad++) {
        if (quadOnAxis[iquad]) continue;
        for (int j = 0; j < 4; j++) {
            int nTag = mConnectivity(iquad, j);
            if (nodeNearAxis[nTag]) mVicinalAxis(iquad, j) = j;
        }
    }
    MultilevelTimer::end("Process Exodus Vicinal", 2);
    
    // check if ocean presents in mesh
    MultilevelTimer::begin("Process Exodus Check Ocean", 2);
    std::string strVs = isIsotropic() ? "VS_0" : "VSV_0";
    for (int iQuad = 0; iQuad < getNumQuads(); iQuad++) {
        if (iQuad % XMPI::nproc() != XMPI::rank()) {
            continue;
        }
        if (getSideSurface(iQuad) == -1) {
            continue;
        }
        double vs = mElementalVariables.at(strVs)(iQuad);
        if (vs < tinyDouble) {
            throw std::runtime_error("ExodusModel::finishReading || "
                "Ocean is detected in mesh. By far, realistic ocean is not implemented in AxiSEM3D. ||"
                "Use non-ocean models in the Mesher and add ocean load in inparam.basic.");
        }
    }
    MultilevelTimer::end("Process Exodus Check Ocean", 2);
    
    // CMB and ICB
    // MultilevelTimer::begin("Process Exodus CMB & ICB", 2);
    // mR_CMB = DBL_MIN;
    // mR_ICB = DBL_MAX;
    // for (int i = 0; i < getNumQuads(); i++) {
    //     bool isFluid = getElementalVariables("fluid", i) > .5;
    //     bool isAxis = getSideAxis(i) >= 0;
    //     if (isFluid && isAxis) {
    //         for (int j = 0; j < 4; j++) {
    //             double s = mNodalS(mConnectivity(i, j));
    //             double z = mNodalZ(mConnectivity(i, j));
    //             double r = std::sqrt(s * s + z * z);
    //             mR_CMB = std::max(mR_CMB, r);
    //             mR_ICB = std::min(mR_ICB, r);
    //         }
    //     }
    // }
    // MultilevelTimer::end("Process Exodus CMB & ICB", 2);
}

std::string ExodusModel::verbose() const {
    std::stringstream ss;
    ss << "\n======================= Exodus Model =======================" << std::endl;
    ss << "  Overview__________________________________________________" << std::endl;
    ss << "    Exodus Title      =   " << mExodusTitle << std::endl;
    ss << "    Mesh CS Type      =   " << (isCartesian() ? "Cartesian" : "Spherical") << std::endl;
    ss << "    Number of Nodes   =   " << getNumNodes() << std::endl;
    ss << "    Number of Quads   =   " << getNumQuads() << std::endl;
    ss << "  Global Variables__________________________________________" << std::endl;
    int widthname = -1;
    for (auto it = mGlobalVariables.begin(); it != mGlobalVariables.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));
    }
    for (auto it = mGlobalRecords.begin(); it != mGlobalRecords.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));    
    }
    for (auto it = mGlobalVariables.begin(); it != mGlobalVariables.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << "   =   " << it->second << std::endl;
    }
    for (auto it = mGlobalRecords.begin(); it != mGlobalRecords.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << "   =   " << it->second << std::endl;
    }
    ss << "  Connectivity______________________________________________" << std::endl;
    int width = (int)std::log10(std::max(getNumQuads(), getNumNodes())) + 1;
    ss << "    " << std::setw(width) << 0 << ": ";
    for (int j = 0; j < 4; j++) {
        ss << std::setw(width) << mConnectivity(0, j) << " ";
    }
    ss << std::endl << "    " << std::setw(width) << "..." << std::endl;
    ss << "    " << std::setw(width) << getNumQuads() - 1 << ": ";
    for (int j = 0; j < 4; j++) {
        ss << std::setw(width) << mConnectivity(getNumQuads() - 1, j) << " ";
    }
    ss << std::endl;
    ss << "  Coordinates_______________________________________________" << std::endl;
    ss << "    " << std::setw(width) << 0 << ": ";
    ss << std::setw(13) << mNodalS(0) << std::setw(13) << mNodalZ(0) << std::endl;
    ss<< "    " << std::setw(width) << "..." << std::endl;
    ss << "    " << std::setw(width) << getNumNodes() - 1 << ": ";
    ss << std::setw(13) << mNodalS(getNumNodes() - 1) << std::setw(13) << mNodalZ(getNumNodes() - 1) << std::endl;
    ss << "  Elemental Variables_______________________________________" << std::endl;
    widthname = -1;
    for (auto it = mElementalVariables.begin(); it != mElementalVariables.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));
    }
    for (auto it = mElementalVariables.begin(); it != mElementalVariables.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << ": ";
        ss << std::setw(13) << it->second(0) << ", ..., ";
        ss << std::setw(13) << it->second(getNumQuads() - 1) << std::endl;        
    }
    ss << "  Side Sets_________________________________________________" << std::endl;
    widthname = -1;
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));
    }
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << ":   ";
        int pair = 0;
        for (int q = 0; q < getNumQuads(); q++) {
            pair += (int)(it->second(q) >= 0);
        }
        ss << pair << " edges" << std::endl;
    }
    ss << "  Miscellaneous_____________________________________________" << std::endl;
    ss << "    Distance Tolerance / m   =   " << mDistTolerance << std::endl;
    if (!isCartesian()) {
        ss << "  External__________________________________________________" << std::endl;
        ss << "    Num. Ellipticity Spline Knots   =   " << mEllipKnots.size() << std::endl;
    }
    ss << "======================= Exodus Model =======================\n" << std::endl;
    return ss.str();
}

void ExodusModel::buildInparam(ExodusModel *&exModel, const Parameters &par, 
    AttParameters *&attPar, int verbose) {
    if (exModel) {
        delete exModel;
    }
    std::string exfile = par.getValue<std::string>("MODEL_1D_EXODUS_MESH_FILE");
    exfile = Parameters::sInputDirectory + "/" + exfile;
    exModel = new ExodusModel(exfile);
    exModel->initialize();
    if (verbose) {
        XMPI::cout << exModel->verbose();
    }
    
    // form attenuation parameters
    if (exModel->hasAttenuation()) {
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
        if (attPar) {
            delete attPar;
        }
        attPar = new AttParameters(nr_lin_solids, f_min, f_max, f_ref, w, y);
    } else {
        if (attPar) {
            delete attPar;
            attPar = 0;
        }
    }
    
    // ellipticity
    if (exModel->isCartesian()) {
        return;
    }
    std::string emode = par.getValue<std::string>("MODEL_3D_ELLIPTICITY_MODE");
    if (boost::iequals(emode, "off")) {
        // no ellipticity
        Geodesy::setup(exModel->getROuter(), 0., RDColX::Zero(0), RDColX::Zero(0));
    } else {
        double inv_f = par.getValue<double>("MODEL_3D_ELLIPTICITY_INVF");
        if (inv_f <= 0.) {
            throw std::runtime_error("ExodusModel::buildInparam || Invalid flattening.");
        } 
        Geodesy::setup(exModel->getROuter(), 1. / inv_f, exModel->mEllipKnots, exModel->mEllipCoeffs);
    }
}

