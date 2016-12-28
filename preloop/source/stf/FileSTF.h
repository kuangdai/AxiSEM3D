// FileSTF.h
// created by Kuangdai on 11-May-2016 
// read stf from file
// file format:
// dt shift series

#pragma once

#include "STF.h"
#include <string>

class FileSTF: public STF {
public: 
    FileSTF(const std::string &fname);
    std::string verbose() const;
    
private:
    std::string mFileName;
};