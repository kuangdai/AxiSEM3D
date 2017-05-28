// WisdomNrField.h
// created by Kuangdai on 10-Oct-2016 
// wisdom nr field

#pragma once
#include "NrField.h"

class NuWisdom;

class WisdomNrField: public NrField {
public:
    WisdomNrField(bool useLucky, const std::string &fname, double factor);
    ~WisdomNrField();
    
    int getNrAtPoint(const RDCol2 &coords) const;
    
    std::string verbose() const;
    
private:
    std::string mFileName;
    double mFactor;
    
    NuWisdom *mNuWisdom;
    const int mNumInterpPoints = 4;
};

