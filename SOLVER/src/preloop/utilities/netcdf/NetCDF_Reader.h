// NetCDF_Reader.h
// created by Kuangdai on 17-May-2017 
// NetCDF Reader

#pragma once

#include <string>
#include "eigenp.h"

class NetCDF_Reader {
public:    
    // file
    virtual void open(const std::string &fname);
    virtual void close();
    bool isOpen() const {return mFileName != "";};
    
    // read
    virtual void readMetaData(const std::string &vname, RDColX &data, std::vector<size_t> &dims) const;
    virtual void read1D(const std::string &vname, RDColX &data) const;
    virtual void read2D(const std::string &vname, RDMatXX &data) const;
    virtual void read3D(const std::string &vname, std::vector<RDMatXX> &data) const;
    
    // string
    void readString(const std::string &vname, std::vector<std::string> &data) const;
    
    // check
    static bool isNetCDF(const std::string &fname);
    
    // determine netcdf or ascii
    static NetCDF_Reader *createOpenNetCDF_Reader(const std::string &fname);
    
private:    
    // error handler
    void netcdfError(const int retval, const std::string &func_name) const;
    
private:
    int mFileID = -1;
    
protected:
    std::string mFileName = "";
};

