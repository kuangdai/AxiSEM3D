// NrFieldEnhance.h
// created by Kuangdai on 2-Feb-2017
// enhanced nr integer field

#pragma once
#include "eigenp.h"

class Parameters;

class NrFieldEnhance {
public:
    virtual ~NrFieldEnhance() {};
    NrFieldEnhance(int ref, bool decrease);
    void updateNrAtPoint(const RDCol2 &sz_target, int nr_base, int &nr_cur) const;
    
    static void buildInparam(std::vector<NrFieldEnhance *> &nrf, 
        const Parameters &par, int verbose);
        
    virtual std::string verbose() const = 0;
    virtual double getValue(const RDCol2 &sz_target) const = 0;
    
protected:
    int mReference;
    bool mDecrease;
};

