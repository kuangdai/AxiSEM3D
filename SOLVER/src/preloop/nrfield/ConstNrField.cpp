// ConstNrField.cpp
// created by Kuangdai on 13-May-2016 
// constant nr integer field

#include "ConstNrField.h"
#include <sstream>

ConstNrField::ConstNrField(bool useLucky, int nu): 
NrField(useLucky), mNu(nu) {
    if (mNu < 0) {
        throw std::runtime_error("ConstNrField::ConstNrField || Negative Nu.");
    }
}

int ConstNrField::getNrAtPoint(const RDCol2 &coords) const {
    int nr = 2 * mNu + 1;
    return nr;
}

std::string ConstNrField::verbose() const {
    std::stringstream ss;
    ss << "\n================= Fourier Expansion Order ==================" << std::endl;
    ss << "  Type                     =   Constant" << std::endl;
    ss << "  Specified Order          =   " << mNu << std::endl;
    ss << "  Use FFTW Lucky Numbers   =   " << (mUseLuckyNumber ? "YES" : "NO") << std::endl;
    ss << "================= Fourier Expansion Order ==================\n" << std::endl;
    return ss.str();
}

