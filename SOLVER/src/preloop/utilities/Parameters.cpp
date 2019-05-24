// Parameters.cpp
// created by Kuangdai on 28-Jun-2016 
// simulation parameters

#include "Parameters.h"
#include "XMPI.h"
#include "global.h"
#include <fstream>
#include <sstream>
#include <iomanip>

std::string Parameters::sInputDirectory;
std::string Parameters::sOutputDirectory;

void Parameters::initReadBcast() {
    registerAll();
    readParFile(Parameters::sInputDirectory + "/inparam.model");
    readParFile(Parameters::sInputDirectory + "/inparam.nu");
    readParFile(Parameters::sInputDirectory + "/inparam.time_src_recv");
    readParFile(Parameters::sInputDirectory + "/inparam.advanced");
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
    registerPar("MODEL_2D_AZIMUTH");
    registerPar("MODEL_PLOT_SLICES_NUM");
    registerPar("MODEL_PLOT_SLICES_LIST");
    registerPar("ATTENUATION");
    
    // inparam.nu
    registerPar("NU_TYPE");
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
    registerPar("NU_WISDOM_LEARN_EPSILON");
    registerPar("NU_WISDOM_LEARN_INTERVAL");
    registerPar("NU_WISDOM_LEARN_OUTPUT");
    registerPar("NU_WISDOM_REUSE_INPUT");
    registerPar("NU_WISDOM_REUSE_FACTOR");
    registerPar("NU_USER_PARAMETER_LIST");
    
    // inparam.time_src_recv
    registerPar("TIME_DELTA_T");
    registerPar("TIME_DELTA_T_FACTOR");
    registerPar("TIME_RECORD_LENGTH");
    registerPar("SOURCE_TYPE");
    registerPar("SOURCE_FILE");
    registerPar("SOURCE_TIME_FUNCTION");
    registerPar("SOURCE_STF_HALF_DURATION");
    registerPar("OUT_STATIONS_FILE");
    registerPar("OUT_STATIONS_DUPLICATED");
    registerPar("OUT_STATIONS_SYSTEM");
    registerPar("OUT_STATIONS_FORMAT");
    registerPar("OUT_STATIONS_COMPONENTS");
    registerPar("OUT_STATIONS_RECORD_INTERVAL");
    registerPar("OUT_STATIONS_DUMP_INTERVAL");
    registerPar("OUT_STATIONS_WHOLE_SURFACE");
    registerPar("OUT_STATIONS_DEPTH_REF");
    
    // inparam.advanced
    registerPar("ATTENUATION_CG4");
    registerPar("ATTENUATION_SPECFEM_LEGACY");
    registerPar("ATTENUATION_QKAPPA");
    registerPar("DD_PROC_INTERVAL");
    registerPar("DD_NCUTS_PER_PROC");
    registerPar("OPTION_VERBOSE_LEVEL");
    registerPar("OPTION_STABILITY_INTERVAL");
    registerPar("OPTION_LOOP_INFO_INTERVAL");
    registerPar("DEVELOP_MAX_TIME_STEPS");
    registerPar("DEVELOP_NON_SOURCE_MODE");
    registerPar("DEVELOP_DIAGNOSE_PRELOOP");
    registerPar("DEVELOP_MEASURED_COSTS");
    registerPar("DEVELOP_RANDOMIZE_DISP0");
    registerPar("FFTW_LUCKY_NUMBER");
    registerPar("FFTW_DISABLE_WISDOM");
    
}

void Parameters::parseLine(const std::string &line) {
    std::vector<std::string> strs = splitString(line, "\t ");
    std::string key = strs[0];
    if (key == "#") {
        return;
    }
    auto it = mKeyValues.find(key);
    if (it != mKeyValues.end()) {
        it->second.insert(it->second.end(), strs.begin() + 1, strs.end());
    }
}

void Parameters::registerPar(const std::string &key) {
    mKeyValues.insert(std::pair<std::string, std::vector<std::string>>
        (key, std::vector<std::string>()));
}

void Parameters::readParFile(const std::string &fname) {
    std::vector<std::string> all_lines;
    if (XMPI::root()) {
        std::string line;
        std::fstream fs(fname, std::fstream::in);
        if (!fs) {
            throw std::runtime_error("Parameters::readParFile || "
                "Error opening parameter file: ||" + fname);
        }
        while (getline(fs, line)) {
            all_lines.push_back(line);
        }
        fs.close();
    }
    XMPI::bcast(all_lines);
    for (int i = 0; i < all_lines.size(); i++) {
        parseLine(all_lines[i]);
    }
}

std::string Parameters::verbose() const {
    std::stringstream ss;
    int width = -1;
    for (auto it = mKeyValues.begin(); it != mKeyValues.end(); it++) {
        width = std::max(width, (int)it->first.length());
    }
        
    ss << "\n======================== Parameters ========================" << std::endl;
    for (auto it = mKeyValues.begin(); it != mKeyValues.end(); it++) {
        ss << "  " << std::setw(width) << std::left << it->first << "   =   ";
        for (const auto &value: it->second) {
            ss << value << " ";
        }
        ss << std::endl;
    }
    // append mpi nproc
    ss << "  " << std::setw(width) << std::left << "mpi nproc" << "   =   ";
    ss << XMPI::nproc() << std::endl;
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
    if (boost::iequals(vstr, "essential")) {
        verbose = 1;
    } else if (boost::iequals(vstr, "detailed")) {
        verbose = 2;
    } else {
        verbose = 0;
    }
    if (verbose) {
        XMPI::cout << par->verbose();
    }
}


std::vector<std::string> Parameters::splitString(const std::string &in, 
    const std::string &sep) {
    // trim
    std::string line(in);
    boost::trim_if(line, boost::is_any_of("\t "));
    // split
    std::vector<std::string> strs;
    boost::split(strs, line, boost::is_any_of(sep), boost::token_compress_on);
    return strs;
}
