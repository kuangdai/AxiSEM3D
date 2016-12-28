// WisdomNrField.cpp
// created by Kuangdai on 10-Oct-2016 
// wisdom nr field

#include "WisdomNrField.h"
#include "XMath.h"
#include <sstream>

#include "NuWisdom.h"

WisdomNrField::WisdomNrField(bool useLucky, const std::string &fname, double factor): 
NrField(useLucky), mFileName(fname), mFactor(factor) {
    mNuWisdom = new NuWisdom();
    mNuWisdom->readFromFile(fname, false);
}

WisdomNrField::~WisdomNrField() {
    delete mNuWisdom;
}

int WisdomNrField::getNrAtPoint(const RDCol2 &coords) const {
    int nu = round(mNuWisdom->getNu(coords(0), coords(1), mNumInterpPoints) * mFactor);
    return nu * 2 + 1;
}

int WisdomNrField::getMaxNr() const {
    int nu = round(mNuWisdom->getMaxNu() * mFactor);
    return nu * 2 + 1;
}

std::string WisdomNrField::verbose() const {
    std::stringstream ss;
    ss << "\n================= Fourier Expansion Order ==================" << std::endl;
    ss << "  Type                     =   Wisdom" << std::endl;
    ss << "  Wisdom File              =   " << mFileName << std::endl;
    ss << "  Wisdom Factor            =   " << mFactor << std::endl;
    ss << "  Compression Ratio        =   " << mNuWisdom->getCompressionRatio() << std::endl;
    ss << "  Maximum Order            =   " << getMaxNr() / 2 << std::endl;
    ss << "  Use FFTW Lucky Numbers   =   " << (mUseLuckyNumber ? "YES" : "NO") << std::endl;
    ss << "================= Fourier Expansion Order ==================\n" << std::endl;
    return ss.str();
}

