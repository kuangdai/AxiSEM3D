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

extern "C" {
    #include "hdf5.h"
};

ExodusModel::ExodusModel(const std::string &fileName): mExodusFileName(fileName) {
    // nothing
}

void ExodusModel::initialize() {
    int io_ws = 0;
    int comp_ws = 8;
    mExodusId = ex_open(mExodusFileName.c_str(), EX_READ, &comp_ws, &io_ws, &mExodusVersion);
    if (mExodusId < 0) throw std::runtime_error("ExodusModel::initialize || "
        "Error opening exodus model file: ||" + mExodusFileName);
    
    int nDim = 0;
    int nBlock = 0;
    int nNodeSets = 0;
    ExodusModel::exodusError(
        ex_get_init(mExodusId, mExodusTitle, &nDim, &mNumNodes, &mNumQuads,
            &nBlock, &nNodeSets, &mNumSideSets), "ex_get_init");   
    
    if (nDim != 2 || nBlock != 1 || nNodeSets != 0)
        throw std::runtime_error("ExodusModel::initialize || Invalid exodus model file: ||" + mExodusFileName);
        
    readGlobalVariables();
    readConnectivity();
    readCoordinates();
    readElementalVariables();
    readSideSets();
    // close as exodus
    ExodusModel::exodusError(ex_close(mExodusId), "ex_close");
    // open as hdf5
    readExternalH5();
    finishReading();
}

void ExodusModel::bcast() {
    XMPI::bcast(mExodusId);
    XMPI::bcast(mExodusVersion);
    XMPI::bcast(mExodusTitle, MAX_LINE_LENGTH + 1);
    XMPI::bcast(mNumNodes);
    XMPI::bcast(mNumQuads);
    XMPI::bcast(mGlobalVariables);
    XMPI::bcast(mGlobalRecords);
    XMPI::bcast(mConnectivity);
    XMPI::bcast(mNodalS);
    XMPI::bcast(mNodalZ);
    XMPI::bcast(mElementalVariables);
    XMPI::bcast(mNumSideSets);
    XMPI::bcast(mSideSets);
    XMPI::bcast(mDistTolerance);
    XMPI::bcast(mROuter);
    XMPI::bcast(mAveGLLSpacing);
    XMPI::bcast(mVicinalAxis);
    XMPI::bcast(mCartesian);
    XMPI::bcast(mSSNameAxis);
    XMPI::bcast(mSSNameSurface);
    XMPI::bcast(mEllipKnots);
    XMPI::bcast(mEllipCoeffs);
}

void ExodusModel::readGlobalVariables() {
    ///// float type /////
    // get number of variables
    int nVars = 0;
    exodusError(ex_get_var_param(mExodusId, "g", &nVars), "ex_get_var_param");
    // variable names 
    char *varNames[nVars];
    for (int i = 0; i < nVars; i++) varNames[i] = new char[256];
    exodusError(ex_get_var_names(mExodusId, "g", nVars, varNames), "ex_get_var_names");
    // variable values
    std::vector<double> varValues(nVars);
    exodusError(ex_get_glob_vars(mExodusId, 1, nVars, varValues.data()), "ex_get_glob_vars");
    // insert to map
    for (int i = 0; i < nVars; i++) {
        std::string gname(varNames[i]);
        // clarify here the meaning of dt
        if (gname == "dt") gname = "dt (nPol = 1)";
        mGlobalVariables.insert(std::pair<std::string, double>(gname, varValues[i]));
        delete [] varNames[i];
    }
    
    ///// string type /////
    // get number of records
    int nRecords = ex_inquire_int(mExodusId, EX_INQ_INFO);
    // read records
    char *records[nRecords];
    for (int i = 0; i < nRecords; i++) records[i] = new char[256];
    exodusError(ex_get_info(mExodusId, records), "ex_get_info");
    // insert to map
    for (int i = 0; i < nRecords; i++) {
        // records to be excluded
        // It may not feel good to see one's name present in the simulation info 
        std::vector<std::string> excluded = {"cmdl 0", "cmdl 1", "host", "python version", "user"};
        std::string recstr(records[i]);
        delete [] records[i];
        std::vector<std::string> substrs;
        boost::trim_if(recstr, boost::is_any_of("\t "));
        boost::split(substrs, recstr, boost::is_any_of("="), boost::token_compress_on);
        if (std::find(excluded.begin(), excluded.end(), boost::trim_copy(substrs[0])) == excluded.end()) {
            mGlobalRecords.insert(std::pair<std::string, std::string>(
                boost::trim_copy(substrs[0]), boost::trim_copy(substrs[1])));
            if (substrs[0].find("crdsys") != std::string::npos) 
                mCartesian = (substrs[1].find("cartesian") != std::string::npos);
        }
    }
    // name of axis and surface sets
    mSSNameAxis = mCartesian ? "x0" : "t0";
    mSSNameSurface = mCartesian ? "y1" : "r1";
}

void ExodusModel::readConnectivity() {
    mConnectivity.resize(mNumQuads);
    exodusError(ex_get_elem_conn(mExodusId, 1, mConnectivity.data()), "ex_get_elem_conn");
    // let start from 0
    for (auto && connect: mConnectivity) {
        connect[0] -= 1;
        connect[1] -= 1;
        connect[2] -= 1;
        connect[3] -= 1;
    }    
}

void ExodusModel::readCoordinates() {
    mNodalS.resize(mNumNodes);
    mNodalZ.resize(mNumNodes);
    exodusError(ex_get_coord(mExodusId, mNodalS.data(), mNodalZ.data(), NULL), "ex_get_coord");
    
    // NOTE: we temporarily treat Cartesian meshes as special cases of spherical meshes
    //       by means of moving it to the "north pole". The introduced global 
    //       curvature should be ignorable, or the problem itself is ill-defined
    //       as a local problem.
    if (mCartesian) {
        double R_EARTH = 6371e3;
        double maxz = -1.;
        for (int i = 0; i < mNumNodes; i++) 
            maxz = std::max(maxz, mNodalZ[i]);
        for (int i = 0; i < mNumNodes; i++) 
            mNodalZ[i] += R_EARTH - maxz;
    }
}


void ExodusModel::readElementalVariables() {
    // get number of variables
    int nVars = 0;
    exodusError(ex_get_var_param(mExodusId, "e", &nVars), "ex_get_var_param");
    // variable names 
    char *varNames[nVars];
    for (int i = 0; i < nVars; i++) varNames[i] = new char[256];
    exodusError(ex_get_var_names(mExodusId, "e", nVars, varNames), "ex_get_var_names");
    // get variables
    std::vector<double> buffer(mNumQuads);
    for (int i = 0; i < nVars; i++) {
        exodusError(ex_get_elem_var(
            mExodusId, 1, (i + 1), 1, mNumQuads, buffer.data()), "ex_get_elem_var");
        mElementalVariables.insert(std::pair<std::string, std::vector<double>>
            (std::string(varNames[i]), buffer));
        delete [] varNames[i];
    }
}

void ExodusModel::readSideSets() {
    // get names of SideSets
    char *ssNames[mNumSideSets];
    for (int i = 0; i < mNumSideSets; i++) ssNames[i] = new char[256];
    exodusError(ex_get_names(mExodusId, EX_SIDE_SET, ssNames), "ex_get_names");
    // get SideSets
    for (int i = 0; i < mNumSideSets; i++) {
        int size, useless;
        exodusError(ex_get_side_set_param(mExodusId, (i + 1), &size, &useless), "ex_get_side_set_param");
        std::vector<int> ebuffer(size);
        std::vector<int> sbuffer(size);
        exodusError(ex_get_side_set(
            mExodusId, (i + 1), ebuffer.data(), sbuffer.data()), "ex_get_elem_var");
        std::vector<int> ss(mNumQuads, -1);
        for (int j = 0; j < size; j++) ss[ebuffer[j] - 1] = sbuffer[j] - 1;
        mSideSets.insert(std::pair<std::string, std::vector<int>>(std::string(ssNames[i]), ss));
        delete [] ssNames[i];
    }
}

void ExodusModel::readExternalH5() {
    // open file
    hid_t file_id = H5Fopen(mExodusFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    ///////// read ellipticity /////////
    if (!mCartesian) {
        // read dimensions
        hid_t dset_id = H5Dopen(file_id, "/ellipticity", H5P_DEFAULT);
        hid_t dspace_id = H5Dget_space(dset_id);
        hsize_t dims[2];
        hdf5Error(H5Sget_simple_extent_dims(dspace_id, dims, NULL), "H5Sget_simple_extent_dims");
        int ellip_dim = dims[1];
        // read meta data
        double *data = new double[2 * ellip_dim];
        hdf5Error(H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "H5Dread");
        hdf5Error(H5Dclose(dset_id), "H5Dclose");
        // copy to memory variables
        for (int i = 0; i < ellip_dim; i++) {
            mEllipKnots.push_back(data[i]);
            mEllipCoeffs.push_back(data[i + ellip_dim]);
        }
        // delete meta data
        delete [] data;
    }
    
    // close file
    hdf5Error(H5Fclose(file_id), "H5Dclose");
}

void ExodusModel::finishReading() {
    // distance tolerance
    mDistTolerance = 1e100;
    for (int i = 0; i < mNumQuads; i++) {
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
        mDistTolerance = std::min({dist0, dist1, dist2, dist3, mDistTolerance});
    }
    
    // surface radius
    mROuter = -1.;
    for (int i = 0; i < mNumNodes; i++) 
        mROuter = std::max(mROuter, mNodalZ[i]);
        
    // average gll spacing
    std::vector<std::vector<int>> refElem(mNumNodes, std::vector<int>());
    for (int i = 0; i < mNumQuads; i++) {
        refElem[mConnectivity[i][0]].push_back(i);
        refElem[mConnectivity[i][1]].push_back(i);
        refElem[mConnectivity[i][2]].push_back(i);
        refElem[mConnectivity[i][3]].push_back(i);
    }
    mAveGLLSpacing = std::vector<double>(mNumNodes, 0.);
    for (int i = 0; i < mNumNodes; i++) {
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
    
    // rotate nodes of axial elements such that side 3 is on axis
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
    
    
    // find elements that are not axial but neighboring axial elements
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

    // check if ocean presents in mesh
    std::string strVs = isIsotropic() ? "VS_0" : "VSV_0";
    for (int iQuad = 0; iQuad < mNumQuads; iQuad++) {
        if (getSideSurface(iQuad) == -1) continue;
        double vs = mElementalVariables.at(strVs)[iQuad];
        if (vs < tinyDouble) throw std::runtime_error("ExodusModel::finishReading || "
            "Ocean is detected in mesh. By far, realistic ocean is not implemented in AxiSEM3D. ||"
            "Use non-ocean models in the Mesher and add ocean load in inparam.basic.");
    }
}

void ExodusModel::exodusError(const int retval, const std::string &func_name) {
    if (retval) throw std::runtime_error("ExodusModel::exodusError || "
        "Error in exodus function: " + func_name);
}

void ExodusModel::hdf5Error(const int retval, const std::string &func_name) {
    if (retval < 0) throw std::runtime_error("ExodusModel::hdf5Error || "
        "Error in exodus function: " + func_name);
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
    std::string exfile = par.getValue<std::string>("MODEL_EXODUS_MESH_FILE");
    exfile = Parameters::sInputDirectory + "/" + exfile;
    exModel = new ExodusModel(exfile);
    if (XMPI::root()) exModel->initialize();
    exModel->bcast();
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
    std::string emode = par.getValue<std::string>("MODEL_ELLIPTICITY_MODE");
    if (boost::iequals(emode, "off")) {
        // no ellipticity
        XMath::setEllipticity(0., exModel->mROuter, std::vector<double>(), std::vector<double>());
    } else {
        double inv_f = par.getValue<double>("MODEL_ELLIPTICITY_INVF");
        if (inv_f <= 0.) throw std::runtime_error("ExodusModel::buildInparam || Invalid flattening."); 
        XMath::setEllipticity(1. / inv_f, exModel->mROuter, exModel->mEllipKnots, exModel->mEllipCoeffs);
    }
}



