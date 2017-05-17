// NetCDF_ReaderAscii.h
// created by Kuangdai on 17-May-2017 
// ascii alternative for NetCDF

#pragma once

#include <fstream>
#include "NetCDF_Reader.h"

class NetCDF_ReaderAscii: public NetCDF_Reader {
public:
    // io
    void open(const std::string &fname);
    void close();
    
    // read
    void readMetaData(const std::string &vname, RDColX &data, std::vector<size_t> &dims);
    void read1D(const std::string &vname, RDColX &data);
    void read2D(const std::string &vname, RDMatXX &data);
    void read3D(const std::string &vname, std::vector<RDMatXX> &data);
    
    // check
    static bool isNetCDFAscii(const std::string &fname);
    
private:
    bool checkVarStart(const std::string &line, const std::string &vname);
    bool checkVarEnd(std::fstream &fs, const std::string &vname);
    void getVarDims(const std::string &line, std::vector<size_t> &dims);
    
private:
    std::fstream mFile;
};

