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
    mPWD = mFileID;
    if (overwrite) {
        netcdfError(nc_enddef(mFileID), "nc_enddef");
    }
}

void NetCDF_Writer::close() {
    if (isOpen()) {
        netcdfError(nc_close(mFileID), "nc_close");
        mPWD = mFileID = -1;
        mFileName = "";
    }
}

void NetCDF_Writer::writeString(const std::string &vname, const std::string &data) const {
    std::vector<size_t> dims;
    dims.push_back(data.length());
    defineVariable(vname, dims, '0');
    int varid = inquireVariable(vname);
    if (data.length() > 0) {
        if (nc_put_var(mPWD, varid, data.c_str()) != NC_NOERR) {
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
        if (nc_put_var(mPWD, varid, buf.data()) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Writer::writeStringInByte || "
                "Error writing variable, variable: " + vname + " || NetCDF file: " + mFileName);
        }
    }
}

void NetCDF_Writer::createGroup(const std::string &gname) const {
    int grpid = -1;
    netcdfError(nc_redef(mFileID), "nc_redef");
    if (nc_def_grp(mPWD, gname.c_str(), &grpid) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::createGroup || "
            "Error defining group: " + gname + " || NetCDF file: " + mFileName);
    }
    netcdfError(nc_enddef(mFileID), "nc_enddef");
}

void NetCDF_Writer::goToGroup(const std::string &gname) {
    int grpid = -1;
    if (nc_inq_grp_ncid(mPWD, gname.c_str(), &grpid) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::goToGroup || "
            "Error finding group: " + gname + " || NetCDF file: " + mFileName);
    }
    mPWD = grpid;
}

void NetCDF_Writer::addAttributeString(const std::string &vname, 
    const std::string &attname, const std::string &attvalue) const {
    int varid = inquireVariable(vname);
    netcdfError(nc_redef(mFileID), "nc_redef");
    if (nc_put_att_text(mPWD, varid, attname.c_str(), attvalue.length(), attvalue.c_str()) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::addAttributeString || "
            "Error adding attribute to variable, variable: " + vname + ", attribute: " + attname  
            + " || NetCDF file: " + mFileName);
    }
    netcdfError(nc_enddef(mFileID), "nc_enddef");
}

int NetCDF_Writer::inquireVariable(const std::string &vname) const {
    int varid = -1;
    if (nc_inq_varid(mPWD, vname.c_str(), &varid) != NC_NOERR) {
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

