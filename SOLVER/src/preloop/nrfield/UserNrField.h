// UserNrField.h
// created by Kuangdai on 11-Jun-2017 
// user-defined nr integer field

#pragma once
#include "NrField.h"
#include <vector>
#include <string>

class UserNrField: public NrField {
public:
    UserNrField(bool useLucky, const std::vector<double> &params);
    
    int getNrAtPoint(const RDCol2 &coords) const;
    
    std::string verbose() const;
    
private:
    
    std::vector<double> mParameters;
};

