// SFCoupling1D.h
// created by Kuangdai on 5-Apr-2016 
// 1D solid-fluid boundary condition

#pragma once

#include "SFCoupling.h"

class SFCoupling1D: public SFCoupling {
public:   
    SFCoupling1D(Real ns, Real nz, Real ns_invmf, Real nz_invmf): 
        mNormalS_unassembled(ns), 
        mNormalZ_unassembled(nz), 
        mNormalS_assembled_invMassFluid(ns_invmf), 
        mNormalZ_assembled_invMassFluid(nz_invmf) {};
    
    // solid-fluid coupling
    void coupleFluidToSolid(const CColX &fluidStiff, CMatX3 &solidStiff) const; 
    void coupleSolidToFluid(const CMatX3 &solidDispl, CColX &fluidStiff) const;
    
    // verbose
    std::string verbose() const {return "SFCoupling1D";};    
    
private:    
    Real mNormalS_unassembled;
    Real mNormalZ_unassembled;
    Real mNormalS_assembled_invMassFluid;
    Real mNormalZ_assembled_invMassFluid;
};
