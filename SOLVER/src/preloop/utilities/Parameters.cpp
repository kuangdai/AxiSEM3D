// Parameters.cpp
// created by Kuangdai on 28-Jun-2016 
// simulation parameters

#include "Parameters.h"
#include "XMPI.h"
#include "global.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

std::string Parameters::sInputDirectory;
std::string Parameters::sOutputDirectory;

void Parameters::initReadBcast() {
    if (XMPI::root()) {
        registerAll();
        readParFile(Parameters::sInputDirectory + "/inparam.model");
        readParFile(Parameters::sInputDirectory + "/inparam.nu");
        readParFile(Parameters::sInputDirectory + "/inparam.time_src_recv");
        readParFile(Parameters::sInputDirectory + "/inparam.advanced");
    }
    XMPI::bcast(mKeyValues);
}

void Parameters::registerAll() {

    // inparam.model
    registerPar("MODEL_1D_EXODUS_MESH_FILE");
    registerPar("MODEL_3D_VOLUMETRIC_NUM");
    registerPar("MODEL_3D_VOLUMETRIC_LIST");
    registerPar("MODEL_3D_GEOMETRIC_NUM");
    registerPar("MODEL_3D_GEOMETRIC_LIST");
    registerPar("MODEL_3D_ELLIPTICITY_MODE");
    registerPar("MODEL_3D_ELLIPTICITY_INVF");
    registerPar("MODEL_3D_OCEAN_LOAD");
    registerPar("MODEL_2D_MODE");
    registerPar("MODEL_2D_LATITUDE");
    registerPar("MODEL_2D_LONGITUDE");
    registerPar("MODEL_PLOT_SLICES_NUM");
    registerPar("MODEL_PLOT_SLICES_LIST");
    registerPar("ATTENUATION");
    
    // inparam.nu
    registerPar("NU_TYPE");
    registerPar("NU_FFTW_LUCKY_NUMBER");
    registerPar("NU_CONST");
    registerPar("NU_EMP_REF");
    registerPar("NU_EMP_MIN");
    registerPar("NU_EMP_SCALE_AXIS");
    registerPar("NU_EMP_POW_AXIS");
    registerPar("NU_EMP_SCALE_THETA");
    registerPar("NU_EMP_POW_THETA");
    registerPar("NU_EMP_FACTOR_PI");
    registerPar("NU_EMP_THETA_START");
    registerPar("NU_EMP_SCALE_DEPTH");
    registerPar("NU_EMP_FACTOR_SURF");
    registerPar("NU_EMP_DEPTH_START");
    registerPar("NU_EMP_DEPTH_END");
    registerPar("NU_WISDOM_LEARN");
    registerPar("NU_WISDOM_LEARN_AIM");
    registerPar("NU_WISDOM_LEARN_INTERVAL");
    registerPar("NU_WISDOM_LEARN_OUTPUT");
    registerPar("NU_WISDOM_REUSE_INPUT");
    registerPar("NU_WISDOM_REUSE_FACTOR");
    
    // inparam.time_src_recv
    registerPar("TIME_DELTA_T");
    registerPar("TIME_RECORD_LENGTH");
    registerPar("SOURCE_TYPE");
    registerPar("SOURCE_FILE");
    registerPar("OUT_STATIONS_FILE");
    registerPar("OUT_STATIONS_SYSTEM");
    registerPar("OUT_STATIONS_FORMAT");
    registerPar("OUT_STATIONS_COMPONENTS");
    registerPar("OUT_STATIONS_RECORD_INTERVAL");
    registerPar("OUT_STATIONS_DUMP_INTERVAL");
    
    // inparam.advanced
    registerPar("ATTENUATION_CG4");
    registerPar("ATTENUATION_SPECFEM_LEGACY");
    registerPar("ATTENUATION_QKAPPA");
    registerPar("DD_BALANCE_ELEMENT_POINT");
    registerPar("DD_NPART_METIS");
    registerPar("DD_COMM_VOL_METIS");
    registerPar("OPTION_VERBOSE_LEVEL");
    registerPar("OPTION_STABILITY_INTERVAL");
    registerPar("OPTION_LOOP_INFO_INTERVAL");
    registerPar("DEVELOP_MAX_TIME_STEPS");
    registerPar("DEVELOP_NON_SOURCE_MODE");
    registerPar("DEVELOP_DIAGNOSE_PRELOOP");
    registerPar("DEVELOP_MEASURED_ELE_COST");
    
}

void Parameters::parseLine(const std::string &line_in) {
    std::string line = line_in;
    std::vector<std::string> strs;
    boost::trim_if(line, boost::is_any_of("\t "));
    boost::split(strs, line, boost::is_any_of("\t "), boost::token_compress_on);
    std::string key = strs[0];
    if (key == "#") return;
    auto it = mKeyValues.find(key);
    if (it != mKeyValues.end()) 
        it->second = std::vector<std::string>(strs.begin() + 1, strs.end());
}

void Parameters::registerPar(const std::string &key) {
    mKeyValues.insert(std::pair<std::string, std::vector<std::string>>
        (key, std::vector<std::string>()));
}

void Parameters::readParFile(const std::string &fname) {
    std::fstream fs(fname, std::fstream::in);
    if (!fs) throw std::runtime_error("Parameters::readParFile || "
        "Error opening parameter file: ||" + fname);
    
    std::string line;
    while (getline(fs, line)) parseLine(line);
}

std::string Parameters::verbose() const {
    std::stringstream ss;
    int width = -1;
    for (auto it = mKeyValues.begin(); it != mKeyValues.end(); it++) 
        width = std::max(width, (int)it->first.length());
        
    ss << "\n======================== Parameters ========================" << std::endl;
    for (auto it = mKeyValues.begin(); it != mKeyValues.end(); it++) {
        ss << "  " << std::setw(width) << std::left << it->first << "   =   ";
        for (const auto &value: it->second) ss << value << " ";
        ss << std::endl;
    }
    ss << "======================== Parameters ========================\n" << std::endl;
    return ss.str();
}

#include <boost/algorithm/string.hpp>
void Parameters::buildInparam(Parameters *&par, int &verbose) {
    if (par) delete par;
    par = new Parameters();
    par->initReadBcast();
    // initialize verbose level
    std::string vstr = par->getValue<std::string>("OPTION_VERBOSE_LEVEL");
    if (boost::iequals(vstr, "essential"))
        verbose = 1;
    else if (boost::iequals(vstr, "detailed"))
        verbose = 2;
    else 
        verbose = 0;    
    if (verbose == 2) XMPI::cout << par->verbose();
}


