// NetCDF_Writer.h
// created by Kuangdai on 17-May-2017 
// NetCDF Writer
// NOTE: this class always overwrites the NetCDF files.

#pragma once

#include <netcdf.h>
#include <vector>
#include <string>
#include <sstream>
#include <memory>
#include <typeinfo>

class NetCDF_Writer {
public:    
    // file
    void open(const std::string &fname, bool overwrite);
    void close();
    bool isOpen() const {return mFileName != "";};
    
    // define
    template<class base_type>
    void defineVariable(const std::string &vname, const std::vector<size_t> &dims, base_type initValue) const {
        // define the dimensions
        netcdfError(nc_redef(mFileID), "nc_redef");
        std::vector<int> dimids;
        size_t total = 1;
        for (int i = 0; i < dims.size(); i++) {
            std::stringstream dimName;
            dimName << "ncdim_" << dims[i];
            int dimid = -1;
            // find ID
            if (nc_inq_dimid(mFileID, dimName.str().c_str(), &dimid) != NC_NOERR) {
                // if not found, create it
                if (nc_def_dim(mFileID, dimName.str().c_str(), dims[i], &dimid) != NC_NOERR) {
                    throw std::runtime_error("NetCDF_Writer::defineVariable || "
                        "Error defining dimension, dimension: " + dimName.str() + " || NetCDF file: " + mFileName);
                }
            }
            dimids.push_back(dimid);
            total *= dims[i];
        }
        
        // define the variable
        int varid = -1;
        nc_type dtype = to_nc_type(initValue);
        if (nc_def_var(mPWD, vname.c_str(), dtype, dims.size(), dimids.data(), &varid) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::defineVariable || "
                "Error defining variable, variable: " + vname + " || NetCDF file: " + mFileName);
        }
        netcdfError(nc_enddef(mFileID), "nc_enddef");
        
        // init value
        std::vector<base_type> constant(total, initValue);
        if (nc_put_var(mPWD, varid, constant.data()) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::defineVariable || "
                "Error initializing variable, variable: " + vname + " || NetCDF file: " + mFileName);
        }
    };
   
    // write an entire variable
    // NOTE: if the Container is an Eigen Matrix, it has to be Eigen::RowMajor
    template<class Container>   
    void writeVariableWhole(const std::string &vname, const Container &data) const {
        int varid = inquireVariable(vname);
        if (data.size() > 0) {
            if (nc_put_var(mPWD, varid, data.data()) != NC_NOERR) {
                throw std::runtime_error("NetCDF_Writer::writeVariableData || "
                    "Error writing variable, variable: " + vname + " || NetCDF file: " + mFileName);
            }
        }
    };
    
    // write a chunk of a variable
    // NOTE: if the Container is an Eigen Matrix, it has to be Eigen::RowMajor
    template<class Container>   
    void writeVariableChunk(const std::string &vname, const Container &data, 
        std::vector<size_t> &start, std::vector<size_t> &count) const {
        int varid = inquireVariable(vname);
        if (data.size() > 0) {
            if (nc_put_vara(mPWD, varid, start.data(), count.data(), data.data()) != NC_NOERR) {
                throw std::runtime_error("NetCDF_Writer::writeVariableChunk || "
                    "Error writing variable, variable: " + vname + " || NetCDF file: " + mFileName);
            }
        }
    };
    
    // string
    void writeString(const std::string &vname, const std::string &data) const;
    void writeStringInByte(const std::string &vname, const std::string &data) const;
    
    // create group
    void createGroup(const std::string &gname) const;
    void goToGroup(const std::string &gname);
    void goToFileRoot() {mPWD = mFileID;};
    
    // add attribute
    void addAttributeString(const std::string &vname, 
        const std::string &attname, const std::string &attvalue) const;
    
    template<class base_type>
    void addAttribute(const std::string &vname, 
        const std::string &attname, base_type attvalue) const {
        int varid = -1;
        int varloc = -1;
        if (vname == "") {
            varid = NC_GLOBAL;
            varloc = mFileID;
        } else {
            varid = inquireVariable(vname);
            varloc = mPWD;
        }
        netcdfError(nc_redef(mFileID), "nc_redef");
        if (nc_put_att(varloc, varid, attname.c_str(), to_nc_type(attvalue), 1, &attvalue) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Reader::addAttribute || "
                "Error adding attribute to variable, variable: " + vname + ", attribute: " + attname  
                + " || NetCDF file: " + mFileName);
        }
        netcdfError(nc_enddef(mFileID), "nc_enddef");    
    };

private:
    // type interpreter
    template<class base_type>
    static nc_type to_nc_type(base_type var) {
        if (typeid(base_type) == typeid(double)) {
            return NC_DOUBLE;
        } else if (typeid(base_type) == typeid(float)) {
            return NC_FLOAT;
        } else if (typeid(base_type) == typeid(int)) {
            return NC_INT;
        } else if (typeid(base_type) == typeid(char)) {
            return NC_CHAR;
        } else if (typeid(base_type) == typeid(signed char)) {
            return NC_BYTE;
        } else if (typeid(base_type) == typeid(long)) {
            return NC_LONG;
        } else if (typeid(base_type) == typeid(long long)) {
            return NC_INT64;
        } else {
            throw std::runtime_error("NetCDF_Writer::to_nc_type || "
                "Error identifying NetCDF type.");
        }
    };
    
    // inquire
    int inquireVariable(const std::string &vname) const;    
    
    // error handler
    void netcdfError(const int retval, const std::string &func_name) const;
    
private:
    int mFileID = -1;
    int mPWD = -1;
    
protected:
    std::string mFileName = "";
};

