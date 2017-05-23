// H5Reader.cpp
// created by Kuangdai on 1-March-2017 
// read hdf5

#include "H5Reader.h"
#include <stdexcept>

// extern "C" {
#include "hdf5.h"
// };

void H5Reader::getAttribute(int fid, const char *key, void *result) {
    hid_t attr = H5Aopen_by_name(fid, ".", key, H5P_DEFAULT, H5P_DEFAULT);
    hid_t atype = H5Aget_type(attr);
    hid_t atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
    hdf5Error(H5Aread(attr, atype_mem, result), "H5Aread");
    hdf5Error(H5Tclose(atype_mem), "H5Tclose");
    hdf5Error(H5Tclose(atype), "H5Tclose");
    hdf5Error(H5Aclose(attr), "H5Tclose");
}

int H5Reader::getDimension1D(int fid, const char *key) {
    hid_t dset_id = H5Dopen(fid, key, H5P_DEFAULT);
    hid_t dspace_id = H5Dget_space(dset_id);
    hsize_t dims[1];
    hdf5Error(H5Sget_simple_extent_dims(dspace_id, dims, NULL), "H5Sget_simple_extent_dims");
    hdf5Error(H5Dclose(dset_id), "H5Dclose");
    return dims[0];
}

void H5Reader::getStringData(int fid, const char *key, int &num, int &len, char *&data, bool alloc) {
    hid_t dset_id = H5Dopen(fid, key, H5P_DEFAULT);
    hid_t dspace_id = H5Dget_space(dset_id);
    hsize_t dims[2];
    dims[0] = dims[1] = 0;
    hdf5Error(H5Sget_simple_extent_dims(dspace_id, dims, NULL), "H5Sget_simple_extent_dims");
    num = dims[0] > 0 ? dims[0] : 1;
    len = dims[1] > 0 ? dims[1] : 1;
    if (alloc && data != 0) {
        delete [] data;
        data = 0;
    }
    if (data == 0) data = new char[num * len + 1];
    data[num * len] = 0;
    hdf5Error(H5Dread(dset_id, H5T_FORTRAN_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "H5Dread");
    hdf5Error(H5Dclose(dset_id), "H5Dclose");
}

void H5Reader::getDoubleData(int fid, const char *key, int &num, int &len, double *&data, bool alloc) {
    hid_t dset_id = H5Dopen(fid, key, H5P_DEFAULT);
    hid_t dspace_id = H5Dget_space(dset_id);
    hsize_t dims[2];
    dims[0] = dims[1] = 0;
    hdf5Error(H5Sget_simple_extent_dims(dspace_id, dims, NULL), "H5Sget_simple_extent_dims");
    num = dims[0] > 0 ? dims[0] : 1;
    len = dims[1] > 0 ? dims[1] : 1;
    // read meta data
    if (alloc && data != 0) {
        delete [] data;
        data = 0;
    }
    if (data == 0) data = new double[num * len];
    hdf5Error(H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "H5Dread");
    hdf5Error(H5Dclose(dset_id), "H5Dclose");
}

void H5Reader::getIntData(int fid, const char *key, int &num, int &len, int *&data, bool alloc) {
    hid_t dset_id = H5Dopen(fid, key, H5P_DEFAULT);
    hid_t dspace_id = H5Dget_space(dset_id);
    hsize_t dims[2];
    dims[0] = dims[1] = 0;
    hdf5Error(H5Sget_simple_extent_dims(dspace_id, dims, NULL), "H5Sget_simple_extent_dims");
    num = dims[0] > 0 ? dims[0] : 1;
    len = dims[1] > 0 ? dims[1] : 1;
    // read meta data
    if (alloc && data != 0) {
        delete [] data;
        data = 0;
    }
    if (data == 0) data = new int[num * len];
    hdf5Error(H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "H5Dread");
    hdf5Error(H5Dclose(dset_id), "H5Dclose");
}

void H5Reader::hdf5Error(const int retval, const std::string &func_name) {
    if (retval < 0) throw std::runtime_error("H5Reader::hdf5Error || "
        "Error in hdf5 function: " + func_name);
}
