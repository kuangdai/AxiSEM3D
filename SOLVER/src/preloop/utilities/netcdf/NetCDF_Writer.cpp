// NetCDF_Writer.h
// created by Kuangdai on 17-May-2017 
// NetCDF Writer
// NOTE: this class always overwrites the NetCDF files.

#include "NetCDF_Writer.h"
#include <sstream>

void NetCDF_Writer::open(const std::string &fname, bool overwrite) {
    close();
    mFileName = fname;
    if (overwrite) {
        if (nc_create(fname.c_str(), NC_NETCDF4, &mFileID) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::open || "
                "Error creating NetCDF file: || " + fname);    
        }
    } else {
        if (nc_open(fname.c_str(),  NC_WRITE | NC_NETCDF4, &mFileID) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::open || "
                "Error opening NetCDF file: || " + fname);    
        }
    }
}

void NetCDF_Writer::close() {
    if (isOpen()) {
        netcdfError(nc_close(mFileID), "nc_close");
        mFileID = -1;
        mFileName = "";
    }
}

void NetCDF_Writer::writeString(const std::string &vname, const std::string &data) const {
    std::vector<size_t> dims;
    dims.push_back(data.length());
    defineVariable(vname, dims, '0');
    int varid = inquireVariable(vname);
    if (data.length() > 0) {
        if (nc_put_var(mFileID, varid, data.c_str()) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::writeString || "
                "Error writing variable, variable: " + vname + " || NetCDF file: " + mFileName);
        }
    }
}

void NetCDF_Writer::writeStringInByte(const std::string &vname, const std::string &data) const {
    std::vector<size_t> dims;
    dims.push_back(data.length());
    defineVariable(vname, dims, (signed char)'0');
    int varid = inquireVariable(vname);
    if (data.length() > 0) {
        std::vector<signed char> buf;
        for (int i = 0; i < data.length(); i++) {
            buf.push_back((signed char)(data.c_str()[i]));
        }
        if (nc_put_var(mFileID, varid, buf.data()) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::writeStringInByte || "
                "Error writing variable, variable: " + vname + " || NetCDF file: " + mFileName);
        }
    }
}

void NetCDF_Writer::createGroup(const std::string &gname) const {
    int grpid = -1;
    nc_redef(mFileID);
    if (nc_def_grp(mFileID, gname.c_str(), &grpid) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::createGroup || "
            "Error defining group: " + gname + " || NetCDF file: " + mFileName);
    }
    nc_enddef(mFileID);
}

int NetCDF_Writer::inquireVariable(const std::string &vname) const {
    int varid = -1;
    if (nc_inq_varid(mFileID, vname.c_str(), &varid) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::inquireVariable || "
            "Error finding variable: " + vname + " || NetCDF file: " + mFileName);
    }
    return varid;
}

void NetCDF_Writer::netcdfError(const int retval, const std::string &func_name) const {
    if (retval != NC_NOERR) {
        throw std::runtime_error("NetCDF_Writer::netcdfError || "
            "Error in NetCDF function: " + func_name + " || NetCDF file: " + mFileName);
    }
}

