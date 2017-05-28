// NrField.h
// created by Kuangdai on 13-May-2016 
// base class of nr integer field

#pragma once
#include "eigenp.h"

class Parameters;

class NrField {
public:
    
    NrField(bool useLucky): mUseLuckyNumber(useLucky) {};
    virtual ~NrField() {};
    
    virtual int getNrAtPoint(const RDCol2 &coords) const = 0;
    
    virtual std::string verbose() const = 0;
    
    static void buildInparam(NrField *&nrf, const Parameters &par, int verbose);
        
    bool useLuckyNumber() const {return mUseLuckyNumber;};
    
protected:
    bool mUseLuckyNumber;
};

