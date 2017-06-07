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
    void readMetaData(const std::string &vname, RDColX &data, std::vector<size_t> &dims) const;
    void read1D(const std::string &vname, RDColX &data) const;
    void read2D(const std::string &vname, RDMatXX &data) const;
    
    // check
    static bool isNetCDFAscii(const std::string &fname);
    
private:
    static bool checkVarStart(const std::string &line, const std::string &vname, const std::string &fname);
    static bool checkVarEnd(std::fstream &fs, const std::string &vname, const std::string &fname);
    static void getVarDims(const std::string &line, std::vector<size_t> &dims);
    
private:
    std::fstream *mFile = 0;
};

