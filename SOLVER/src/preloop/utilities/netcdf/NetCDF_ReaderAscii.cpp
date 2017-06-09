// NetCDF_ReaderAscii.cpp
// created by Kuangdai on 17-May-2017 
// ascii alternative for NetCDF

#include "NetCDF_ReaderAscii.h"
#include <limits>
#include <boost/algorithm/string.hpp>

void NetCDF_ReaderAscii::open(const std::string &fname) {
    close();
    mFileName = fname;
    mFile = new std::fstream(mFileName, std::fstream::in);
    if (!(*mFile)) {
        throw std::runtime_error("NetCDF_ReaderAscii::open || "
            "Error opening NetCDF-alternative ascii file: || " + fname);
    }
}

void NetCDF_ReaderAscii::close() {
    if (isOpen()) {
        (*mFile).close();
        delete mFile;
        mFile = 0;
        mFileName = "";
    }
}

bool NetCDF_ReaderAscii::isNetCDFAscii(const std::string &fname) {
    std::fstream fs(fname, std::fstream::in);
    if (!fs) {
        return false;
    }
    std::string line;
    if (!getline(fs, line)) {
        fs.close();
        return false;
    }
    fs.close();
    return boost::iequals(boost::trim_copy(line), "@@ netcdf ascii format @@");
}

bool NetCDF_ReaderAscii::checkVarStart(const std::string &line, const std::string &vname, 
    const std::string &fname) {
    std::stringstream ss(line);
    
    // find @VAR flag
    std::string flag_target = "@VAR_START";
    std::string flag;
    if (ss >> flag) {
        if (flag != flag_target) {
            return false;
        }
    } else {
        return false;
    }
    
    // find var name
    std::string varName;
    if (ss >> varName) {
        if (varName != vname) {
            return false;
        }
    } else {
        throw std::runtime_error("NetCDF_ReaderAscii::checkVarStart || "
            "No variable name provided after " + flag_target
            + " || NetCDF-alternative ascii file: " + fname);
    }
    
    return true;
}

bool NetCDF_ReaderAscii::checkVarEnd(std::fstream &fs, const std::string &vname, 
    const std::string &fname) {
    // find @VAR flag
    std::string flag_target = "@VAR_END";
    std::string flag;
    if (fs >> flag) {
        if (flag != flag_target) {
            return false;
        }
    } else {
        return false;
    }
    
    // find var name
    std::string varName;
    if (fs >> varName) {
        if (varName != vname) {
            return false;
        }
    } else {
        throw std::runtime_error("NetCDF_ReaderAscii::checkVarEnd || "
            "No variable name provided after " + flag_target
            + " || NetCDF-alternative ascii file: " + fname);
    }
    
    fs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return true;
}

void NetCDF_ReaderAscii::getVarDims(const std::string &line, std::vector<size_t> &dims) {
    std::stringstream ss(line);
    dims.clear();
    size_t dim;
    while (ss >> dim) {
        dims.push_back(dim);
    }
}
