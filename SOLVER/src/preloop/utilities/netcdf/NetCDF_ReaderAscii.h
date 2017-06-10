// NetCDF_ReaderAscii.h
// created by Kuangdai on 17-May-2017 
// ascii alternative for NetCDF

#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

class NetCDF_ReaderAscii {
public:
    // io
    void open(const std::string &fname);
    void close();
    bool isOpen() const {return mFileName != "";};
    
    // read
    template<class Container>
    void readMetaData(const std::string &vname, Container &data, std::vector<size_t> &dims) const {
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
            data.resize(total_len);
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
    };
    
    template<class Container>
    void read1D(const std::string &vname, Container &data) const {
        // read meta data
        std::vector<size_t> dims;
        readMetaData(vname, data, dims);
        
        // check ndims
        int var_ndims = dims.size();
        if (!(var_ndims == 1 || (var_ndims == 2 && dims[0] == 1))) {
            throw std::runtime_error("NetCDF_ReaderAscii::read1D || "
                "Variable is not 1D, Variable = " + vname + " || NetCDF file: " + mFileName);
        }
    };
    
    template<class Container>
    void read2D(const std::string &vname, Container &data) const {
        // read meta data
        std::vector<size_t> dims;
        std::vector<typename Container::Scalar> mdata;
        readMetaData(vname, mdata, dims);
        
        // check ndims
        int var_ndims = dims.size();
        if (var_ndims != 2) {
            throw std::runtime_error("NetCDF_ReaderAscii::read2D || "
                "Variable is not 2D, Variable = " + vname + " || NetCDF file: " + mFileName);
        }
        
        // Container can be both RowMajor and ColMajor
        int pos = 0;
        data = Container::Zero(dims[0], dims[1]);
        for (int j = 0; j < dims[0]; j++) {
            for (int k = 0; k < dims[1]; k++) {
                data(j, k) = mdata[pos++];
            }
        }
    };
    
    void readString(const std::string &vname, std::vector<std::string> &data) const {
        throw std::runtime_error("NetCDF_ReaderAscii::readString || "
            "Not implemented for strings.");
    }
    
    // check
    static bool isNetCDFAscii(const std::string &fname);
    
private:
    static bool checkVarStart(const std::string &line, const std::string &vname, const std::string &fname);
    static bool checkVarEnd(std::fstream &fs, const std::string &vname, const std::string &fname);
    static void getVarDims(const std::string &line, std::vector<size_t> &dims);
    
private:
    std::fstream *mFile = 0;
    std::string mFileName = "";
};

