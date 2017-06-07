// NetCDF_ReaderAscii.cpp
// created by Kuangdai on 17-May-2017 
// ascii alternative for NetCDF

#include "NetCDF_ReaderAscii.h"
#include <stdexcept>
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

void NetCDF_ReaderAscii::readMetaData(const std::string &vname, RDColX &data, std::vector<size_t> &dims) const {
    // rewind
    (*mFile).seekg(0, (*mFile).beg);
    std::string line;
    while (getline((*mFile), line)) {
        // start
        if (!checkVarStart(line, vname, mFileName)) {
            continue;
        }
        
        // dim
        if (!getline((*mFile), line)) {
            throw std::runtime_error("NetCDF_ReaderAscii::readMetaData || "
                "No dimension line after @VAR_START line, varaible: " + vname
                + " || NetCDF-alternative ascii file: " + mFileName);
        }
        getVarDims(line, dims);
        size_t total_len = 1;
        for (int i = 0; i < dims.size(); i++) {
            total_len *= dims[i];
        }
        
        // data
        data = RDColX::Zero(total_len);
        int pos = 0;
        while (pos < total_len) {
            if (! ((*mFile) >> data[pos++])) {
                throw std::runtime_error("NetCDF_ReaderAscii::readMetaData || "
                    "Insufficient data or invalid number format, Variable = " + vname
                    + " || NetCDF-alternative ascii file: " + mFileName);
            }
        }
        
        // end
        if (!checkVarEnd((*mFile), vname, mFileName)) {
            throw std::runtime_error("NetCDF_ReaderAscii::readMetaData || "
                "Error detecting @VAR_END line, Variable = " + vname
                + " || Too many data or invalid number format or missing @VAR_END"
                + " || NetCDF-alternative ascii file: " + mFileName);
        }
        
        return;
    }
    
    throw std::runtime_error("NetCDF_ReaderAscii::readMetaData || "
        "Error detecting @VAR_START line, || "
        "Variable not found, Variable = " + vname
        + " || NetCDF-alternative ascii file: " + mFileName);
}

void NetCDF_ReaderAscii::read1D(const std::string &vname, RDColX &data) const {
    // read meta data
    std::vector<size_t> dims;
    RDColX mdata;
    readMetaData(vname, mdata, dims);
    
    // check ndims
    int var_ndims = dims.size();
    if (var_ndims != 1) {
        throw std::runtime_error("NetCDF_ReaderAscii::read1D || "
            "Variable is not 1D, Variable = " + vname 
            + " || NetCDF-alternative ascii file: " + mFileName);
    }
    
    // get data
    data = mdata;
}

void NetCDF_ReaderAscii::read2D(const std::string &vname, RDMatXX &data) const {
    // read meta data
    std::vector<size_t> dims;
    RDColX mdata;
    readMetaData(vname, mdata, dims);
    
    // check ndims
    int var_ndims = dims.size();
    if (var_ndims != 2) {
        throw std::runtime_error("NetCDF_ReaderAscii::read2D || "
            "Variable is not 2D, Variable = " + vname 
            + " || NetCDF-alternative ascii file: " + mFileName);
    }
    
    // get data
    int pos = 0;
    data = RDMatXX::Zero(dims[0], dims[1]);
    for (int j = 0; j < dims[0]; j++) {
        for (int k = 0; k < dims[1]; k++) {
            data(j, k) = mdata(pos++);
        }
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
