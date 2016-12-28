// NrField.h
// created by Kuangdai on 13-May-2016 
// base class of nr integer field

#pragma once
#include "eigenp.h"

class Parameters;

class NrField {
public:
    virtual ~NrField() {};
    
    NrField(bool useLucky): mUseLuckyNumber(useLucky) {};
    
    virtual int getNrAtPoint(const RDCol2 &coords) const = 0;
    
    virtual int getMaxNr() const = 0;
    
    virtual std::string verbose() const = 0;
    
    static void buildInparam(NrField *&nrf, const Parameters &par, 
        double router, int verbose);
        
    virtual bool useLuckyNumber() const {return mUseLuckyNumber;};
    
protected:
    bool mUseLuckyNumber;
};

