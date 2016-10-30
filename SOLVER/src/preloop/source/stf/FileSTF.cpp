// FileSTF.cpp
// created by Kuangdai on 11-May-2016 
// read stf from file
// file format:
// dt shift series

#include "FileSTF.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

FileSTF::FileSTF(const std::string &fname): 
mFileName(fname) {
    std::fstream fs(fname, std::fstream::in);
    if (!fs) throw std::runtime_error("FileSTF::FileSTF || "
        "Error opening STF data file: ||" + fname);
    fs >> mDeltaT;
    fs >> mShift;
    double x;
    while (fs >> x) mSTF.push_back(x);
    fs.close();
}

std::string FileSTF::verbose() const {
    std::stringstream ss;
    ss << "\n=================== Source Time Function ===================" << std::endl;
    ss << "  Time Step               =   " << mDeltaT << std::endl;
    ss << "  Number of Steps         =   " << mSTF.size() << std::endl;
    ss << "  Total Duration          =   " << mDeltaT * mSTF.size() << std::endl;
    ss << "  Duration after Origin   =   " << mDeltaT * mSTF.size() - mShift << std::endl;
    ss << "  Shift before Origin     =   " << mShift << std::endl;
    ss << "  Time Series Type        =   File based" << std::endl;
    ss << "  File Name               =   " << mFileName << std::endl;
    ss << "=================== Source Time Function ===================\n" << std::endl;
    return ss.str();
}