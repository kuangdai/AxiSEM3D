// NetCDF_Writer.h
// created by Kuangdai on 17-May-2017 
// NetCDF Writer
// NOTE: this class always overwrites the NetCDF files.

#pragma once

#include <string>
#include "eigenp.h"

class NetCDF_Writer {
public:    
    // file
    void open(const std::string &fname);
    void close();
    bool isOpen() const {return mFileName != "";};
    
    // write
    void writeMetaData(const std::string &vname, const RDColX &data, 
        const std::vector<size_t> &dims);
    void write1D(const std::string &vname, const RDColX &data);
    void write2D(const std::string &vname, const RDMatXX &data);
    void write3D(const std::string &vname, const std::vector<RDMatXX> &data);
    
private:    
    // error handler
    void netcdfError(const int retval, const std::string &func_name) const;
    
private:
    int mFileID = -1;
    
protected:
    std::string mFileName = "";
};

