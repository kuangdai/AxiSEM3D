// ConstNrField.h
// created by Kuangdai on 13-May-2016 
// constant nr integer field

#pragma once
#include "NrField.h"

class ConstNrField: public NrField {
public:
    ConstNrField(bool useLucky, int nu);
    
    int getNrAtPoint(const RDCol2 &coords) const;
    
    std::string verbose() const;
    
private:
    
    int mNu;
};

