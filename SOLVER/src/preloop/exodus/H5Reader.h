// H5Reader.h
// created by Kuangdai on 1-March-2017 
// read hdf5

#pragma once

#include <string>

class H5Reader {
public:
    static void getAttribute(int fid, const char *key, void *result);
    static int getDimension1D(int fid, const char *key);
    static void getStringData(int fid, const char *key, int &num, int &len, char *&data, bool alloc = true);
    static void getDoubleData(int fid, const char *key, int &row, int &col, double *&data, bool alloc = true);
    static void getIntData(int fid, const char *key, int &row, int &col, int *&data, bool alloc = true);
    static void hdf5Error(const int retval, const std::string &func_name);
};

