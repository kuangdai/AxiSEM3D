// NetCDF_Writer.h
// created by Kuangdai on 17-May-2017 
// NetCDF Writer
// NOTE: this class always overwrites the NetCDF files.

#include "NetCDF_Writer.h"
#include <stdexcept>
#include <netcdf.h>
#include <sstream>

void NetCDF_Writer::open(const std::string &fname) {
    close();
    mFileName = fname;
    if (nc_create(fname.c_str(), NC_CLOBBER, &mFileID) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Writer::open || "
            "Error creating NetCDF file: || " + fname);    
    }
}

void NetCDF_Writer::close() {
    if (isOpen()) {
        netcdfError(nc_close(mFileID), "nc_close");
        mFileID = -1;
        mFileName = "";
    }
}

void NetCDF_Writer::writeMetaData(const std::string &vname, const RDColX &data, 
    const std::vector<size_t> &dims) const {
    // define the dimensions
    nc_redef(mFileID);
    std::vector<int> dimids;
    for (int i = 0; i < dims.size(); i++) {
        std::stringstream strname;
        strname << "ncdim_" << vname << "_" << i;
        int dimid;
        if (nc_def_dim(mFileID, strname.str().c_str(), dims[i], &dimid)) {
            throw std::runtime_error("NetCDF_Writer::writeMetaData || "
                "Error defining dimensions, variable: " + vname + " || NetCDF file: " + mFileName);
        }
        dimids.push_back(dimid);
    }
    
   // define the variable
   int varid;
   if (nc_def_var(mFileID, vname.c_str(), NC_DOUBLE, dims.size(), dimids.data(), &varid)) {
       throw std::runtime_error("NetCDF_Writer::writeMetaData || "
           "Error defining variable, variable: " + vname + " || NetCDF file: " + mFileName);
   }
   nc_enddef(mFileID);
   
   // put variable data
   if (nc_put_var_double(mFileID, varid, &data[0])) {
       throw std::runtime_error("NetCDF_Writer::writeMetaData || "
           "Error writing variable, variable: " + vname + " || NetCDF file: " + mFileName);
   }
}

void NetCDF_Writer::write1D(const std::string &vname, const RDColX &data) const {
    std::vector<size_t> dims;
    dims.push_back(data.rows());
    writeMetaData(vname, data, dims);
}

void NetCDF_Writer::write2D(const std::string &vname, const RDMatXX &data) const {
    std::vector<size_t> dims;
    dims.push_back(data.rows());
    dims.push_back(data.cols());
    RDColX mdata(dims[0] * dims[1]);
    int pos = 0;
    for (int j = 0; j < dims[0]; j++) {
        for (int k = 0; k < dims[1]; k++) {
            mdata(pos++) = data(j, k);
        }
    }
    writeMetaData(vname, mdata, dims);
}

void NetCDF_Writer::write3D(const std::string &vname, const std::vector<RDMatXX> &data) const {
    std::vector<size_t> dims;
    dims.push_back(data.size());
    dims.push_back(data[0].rows());
    dims.push_back(data[0].cols());
    RDColX mdata(dims[0] * dims[1] * dims[2]);
    int pos = 0;
    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            for (int k = 0; k < dims[2]; k++) {
                mdata(pos++) = data[i](j, k);
            }
        }
    }
    writeMetaData(vname, mdata, dims);
}

void NetCDF_Writer::netcdfError(const int retval, const std::string &func_name) const {
    if (retval != NC_NOERR) throw std::runtime_error("NetCDF_Writer::netcdfError || "
        "Error in NetCDF function: " + func_name + " || NetCDF file: " + mFileName);
}

