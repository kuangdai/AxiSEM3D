// NetCDF_Reader.h
// created by Kuangdai on 17-May-2017 
// NetCDF Reader

#include "NetCDF_Reader.h"
#include <stdexcept>
#include <netcdf.h>

void NetCDF_Reader::open(const std::string &fname) {
    close();
    mFileName = fname;
    if (nc_open(fname.c_str(), 0, &mFileID) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::open || "
            "Error opening NetCDF file: || " + fname);
    }
}

void NetCDF_Reader::close() {
    if (isOpen()) {
        netcdfError(nc_close(mFileID), "nc_close");
        mFileID = -1;
        mFileName = "";
    }
}

void NetCDF_Reader::readMetaData(const std::string &vname, RDColX &data, std::vector<size_t> &dims) const {
    // access variable
    int var_id;
    if (nc_inq_varid(mFileID, vname.c_str(), &var_id) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::readMetaData || "
            "Error finding variable: " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get ndims
    int var_ndims;
    netcdfError(nc_inq_varndims(mFileID, var_id, &var_ndims), "nc_inq_varndims");
    
    // get dim length
    std::vector<int> var_dimids(var_ndims, 0);
    dims.resize(var_ndims);
    netcdfError(nc_inq_vardimid(mFileID, var_id, var_dimids.data()), "nc_inq_vardimid");
    size_t total_len = 1;
    for (int i = 0; i < var_ndims; i++) {
        netcdfError(nc_inq_dimlen(mFileID, var_dimids[i], &dims[i]), "nc_inq_dimlen");
        total_len *= dims[i];
    }
    
    // get data
    data = RDColX::Zero(total_len);
    netcdfError(nc_get_var_double(mFileID, var_id, data.data()), "nc_get_var_double");
}

void NetCDF_Reader::read1D(const std::string &vname, RDColX &data) const {
    // read meta data
    std::vector<size_t> dims;
    RDColX mdata;
    readMetaData(vname, mdata, dims);
    
    // check ndims
    int var_ndims = dims.size();
    if (var_ndims != 1) {
        throw std::runtime_error("NetCDF_Reader::read1D || "
            "Variable is not 1D, Variable = " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get data
    data = mdata;
}

void NetCDF_Reader::read2D(const std::string &vname, RDMatXX &data) const {
    // read meta data
    std::vector<size_t> dims;
    RDColX mdata;
    readMetaData(vname, mdata, dims);
    
    // check ndims
    int var_ndims = dims.size();
    if (var_ndims != 2) {
        throw std::runtime_error("NetCDF_Reader::read2D || "
            "Variable is not 2D, Variable = " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get data
    int pos = 0;
    data = RDMatXX::Zero(dims[0], dims[1]);
    for (int j = 0; j < dims[0]; j++) {
        for (int k = 0; k < dims[1]; k++) {
            data(j, k) = mdata(pos++);
        }
    }
}

void NetCDF_Reader::read3D(const std::string &vname, std::vector<RDMatXX> &data) const {
    // read meta data
    std::vector<size_t> dims;
    RDColX mdata;
    readMetaData(vname, mdata, dims);
    
    // check ndims
    int var_ndims = dims.size();
    if (var_ndims != 3) {
        throw std::runtime_error("NetCDF_Reader::read3D || "
            "Variable is not 3D, Variable = " + vname + " || NetCDF file: " + mFileName);
    }
    
    // structured
    data.clear();
    int pos = 0;
    for (int i = 0; i < dims[0]; i++) {
        RDMatXX mat(dims[1], dims[2]);
        for (int j = 0; j < dims[1]; j++) {
            for (int k = 0; k < dims[2]; k++) {
                mat(j, k) = mdata(pos++);
            }
        }
        data.push_back(mat);
    }
}

void NetCDF_Reader::readMetaData(const std::string &vname, IColX &data, std::vector<size_t> &dims) const {
    // access variable
    int var_id;
    if (nc_inq_varid(mFileID, vname.c_str(), &var_id) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::readMetaData || "
            "Error finding variable: " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get ndims
    int var_ndims;
    netcdfError(nc_inq_varndims(mFileID, var_id, &var_ndims), "nc_inq_varndims");
    
    // get dim length
    std::vector<int> var_dimids(var_ndims, 0);
    dims.resize(var_ndims);
    netcdfError(nc_inq_vardimid(mFileID, var_id, var_dimids.data()), "nc_inq_vardimid");
    size_t total_len = 1;
    for (int i = 0; i < var_ndims; i++) {
        netcdfError(nc_inq_dimlen(mFileID, var_dimids[i], &dims[i]), "nc_inq_dimlen");
        total_len *= dims[i];
    }
    
    // get data
    data = IColX::Zero(total_len);
    netcdfError(nc_get_var_int(mFileID, var_id, data.data()), "nc_get_var_int");
}

void NetCDF_Reader::read1D(const std::string &vname, IColX &data) const {
    // read meta data
    std::vector<size_t> dims;
    IColX mdata;
    readMetaData(vname, mdata, dims);
    
    // check ndims
    int var_ndims = dims.size();
    if (var_ndims != 1) {
        throw std::runtime_error("NetCDF_Reader::read1D || "
            "Variable is not 1D, Variable = " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get data
    data = mdata;
}

void NetCDF_Reader::read2D(const std::string &vname, IMatXX &data) const {
    // read meta data
    std::vector<size_t> dims;
    IColX mdata;
    readMetaData(vname, mdata, dims);
    
    // check ndims
    int var_ndims = dims.size();
    if (var_ndims != 2) {
        throw std::runtime_error("NetCDF_Reader::read2D || "
            "Variable is not 2D, Variable = " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get data
    int pos = 0;
    data = IMatXX::Zero(dims[0], dims[1]);
    for (int j = 0; j < dims[0]; j++) {
        for (int k = 0; k < dims[1]; k++) {
            data(j, k) = mdata(pos++);
        }
    }
}

void NetCDF_Reader::readString(const std::string &vname, std::vector<std::string> &data) const {
    // access variable
    int var_id;
    if (nc_inq_varid(mFileID, vname.c_str(), &var_id) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::readString || "
            "Error finding variable: " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get ndims
    int var_ndims;
    netcdfError(nc_inq_varndims(mFileID, var_id, &var_ndims), "nc_inq_varndims");
    if (var_ndims != 2) {
        throw std::runtime_error("NetCDF_Reader::readString || "
            "Number of dimensions is not 2, Variable = " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get dim length
    std::vector<int> var_dimids(var_ndims, 0);
    std::vector<size_t> dims = std::vector<size_t>(var_ndims, 0);
    netcdfError(nc_inq_vardimid(mFileID, var_id, var_dimids.data()), "nc_inq_vardimid");
    for (int i = 0; i < var_ndims; i++) {
        netcdfError(nc_inq_dimlen(mFileID, var_dimids[i], &dims[i]), "nc_inq_dimlen");
    }
    
    // size
    int numString = dims[0];
    int lenString = dims[1];
    char *cstr = new char[numString * lenString];
    
    netcdfError(nc_get_var_text(mFileID, var_id, cstr), "nc_get_var_text");
    for (int i = 0; i < numString; i++) {
        data.push_back(std::string(&cstr[i * lenString]));
    }
    
    delete [] cstr;
}

bool NetCDF_Reader::isNetCDF(const std::string &fname) {
    int fid, stat;
    stat = nc_open(fname.c_str(), 0, &fid);
    nc_close(fid);
    return stat == NC_NOERR;
}

void NetCDF_Reader::netcdfError(const int retval, const std::string &func_name) const {
    if (retval != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::netcdfError || "
            "Error in NetCDF function: " + func_name + " || NetCDF file: " + mFileName);
    }
}

#include "NetCDF_ReaderAscii.h"
NetCDF_Reader *NetCDF_Reader::createOpenNetCDF_Reader(const std::string &fname) {
    if (isNetCDF(fname)) {
        NetCDF_Reader *reader = new NetCDF_Reader();
        reader->open(fname);
        return reader; 
    }
    
    if (NetCDF_ReaderAscii::isNetCDFAscii(fname)) {
        NetCDF_Reader *reader = new NetCDF_ReaderAscii();
        reader->open(fname);
        return reader; 
    }
    
    throw std::runtime_error("NetCDF_ReaderAscii::createOpenNetCDF_Reader || "
        "Error opening NetCDF or NetCDF-alternative ascii file: || " + fname);
}
