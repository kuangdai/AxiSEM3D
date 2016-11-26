// SlicePlot.cpp
// created by Kuangdai on 26-Nov-2016 
// plot a slice

#include "SlicePlot.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <sstream>
#include "eigenp.h"
#include "XMath.h"
#include "XMPI.h"
#include "Mesh.h"
#include "Parameters.h"

SlicePlot::SlicePlot(const std::string &params, const Mesh *mesh):
mMesh(mesh) {
    std::string mstr = params;
    std::vector<std::string> strs;
    boost::trim_if(mstr, boost::is_any_of("\t "));
    boost::split(strs, mstr, boost::is_any_of("$"), boost::token_compress_on);
    mParName = strs[0];
    bool onSlice = isSlicePar(mParName);
    try {
        if (onSlice) {
            if (strs.size() < 5) throw std::runtime_error("");
            // phi
            mLat = boost::lexical_cast<double>(strs[1]);
            mLon = boost::lexical_cast<double>(strs[2]);
            RDCol3 rtpG;
            rtpG(0) = 1.;
            rtpG(1) = XMath::lat2Theta(mLat, 0.);
            rtpG(2) = XMath::lon2Phi(mLon);
            const RDCol3 &rtpS = XMath::rotateGlob2Src(rtpG, mMesh->mSrcLat, mMesh->mSrcLon, mMesh->mSrcDep);
            mPhi = rtpS(2);
            
            // property type
            if (boost::iequals(strs[3], "1D")) {
                mPropertyType = PropertyTypes::Property1D;
            } else if (boost::iequals(strs[3], "3D")) {
                mPropertyType = PropertyTypes::Property3D;
            } else if (boost::iequals(strs[3], "perturb")) {
                mPropertyType = PropertyTypes::PropertyPerturb;
            } else {
                throw std::runtime_error("");
            }
            
            // sample type
            if (boost::iequals(strs[4], "center")) {
                mSampleType = SampleTypes::Center;
            } else if (boost::iequals(strs[4], "vertex")) {
                mSampleType = SampleTypes::Vertex;
            } else if (boost::iequals(strs[4], "gllpnt")) {
                mSampleType = SampleTypes::GLLPnt;
            } else {
                throw std::runtime_error("");
            }
            
            // radial range
            mRMin = 0.;
            mRMax = 1e100;
            if (strs.size() >= 6) mRMin = boost::lexical_cast<double>(strs[5]);
            if (strs.size() >= 7) mRMin = boost::lexical_cast<double>(strs[6]);
        } else {
            // radial range
            mRMin = 0.;
            mRMax = 1e100;
            if (strs.size() >= 2) mRMin = boost::lexical_cast<double>(strs[1]);
            if (strs.size() >= 3) mRMin = boost::lexical_cast<double>(strs[2]);
        }
        mRMin *= 1e3;
        mRMax *= 1e3;
    } catch (std::exception) {
        throw std::runtime_error("SlicePlot::SlicePlot || "
            "Invalid parameter following ParSeries " + mParName + ".");
    }  
}

std::string SlicePlot::verbose() const {
    std::stringstream ss;
    ss << mParName << " ";
    if (isSlicePar(mParName)) {
        ss << mLat << " " << mLon << " ";
        if (mPropertyType == PropertyTypes::Property1D) ss << "1D ";
        if (mPropertyType == PropertyTypes::Property3D) ss << "3D ";
        if (mPropertyType == PropertyTypes::PropertyPerturb) ss << "perturb ";
        if (mSampleType == SampleTypes::Center) ss << "center ";
        if (mSampleType == SampleTypes::Vertex) ss << "vertex ";
        if (mSampleType == SampleTypes::GLLPnt) ss << "gLLPnt ";
    } else {
        ss << "n.a. n.a. n.a. n.a. ";
    }
    ss << mRMin / 1e3 << " " << mRMax / 1e3;
    return ss.str();
}

bool SlicePlot::isSlicePar(const std::string &parName) const {
    if (boost::iequals(parName, "vp")) return true;
    if (boost::iequals(parName, "vpv")) return true;
    if (boost::iequals(parName, "vph")) return true;
    if (boost::iequals(parName, "vs")) return true;
    if (boost::iequals(parName, "vsv")) return true;
    if (boost::iequals(parName, "vsh")) return true;
    if (boost::iequals(parName, "rho")) return true;
    if (boost::iequals(parName, "eta")) return true;
    if (boost::iequals(parName, "Qmu")) return true;
    if (boost::iequals(parName, "Qkappa")) return true;
    if (boost::iequals(parName, "undulation")) return true;
    
    if (boost::iequals(parName, "nu")) return false;
    if (boost::iequals(parName, "eleType")) return false;
    if (boost::iequals(parName, "eleCost")) return false;
    if (boost::iequals(parName, "unweighted")) return false;
    if (boost::iequals(parName, "weighted")) return false;
    
    throw std::runtime_error("SlicePlot::isSlicePar || Unknown parameter name: " + parName + ".");
}

void SlicePlot::buildInparam(std::vector<SlicePlot *> &splots, 
    const Parameters &par, const Mesh *mesh) {
    // clear
    for (const auto &sp: splots) delete sp;    
    splots.clear();
    // check size
    int nplots = par.getValue<int>("MODEL_PLOT_SLICES_NUM");
    int nsize = par.getSize("MODEL_PLOT_SLICES_LIST");
    if (nplots > nsize) throw std::runtime_error("SlicePlot::buildInparam || "
        "Not enough ParSeries provided in MODEL_PLOT_SLICES_LIST ||"
        "MODEL_PLOT_SLICES_NUM = " + par.getValue<std::string>("MODEL_PLOT_SLICES_NUM") + ".");
    // add to vector
    for (int i = 0; i < nplots; i++) {
        SlicePlot *sp = new SlicePlot(par.getValue<std::string>("MODEL_PLOT_SLICES_LIST", i), mesh);
        splots.push_back(sp);
    }
    
    // verbose
    int verbose = 0;
    std::string vstr = par.getValue<std::string>("OPTION_VERBOSE_LEVEL");
    if (boost::iequals(vstr, "essential")) verbose = 1;
    if (boost::iequals(vstr, "detailed")) verbose = 2;
    if (verbose && nplots > 0) {
        std::stringstream ss;
        ss << "\n======================== Slice Plots =======================" << XMPI::endl;
        ss << "  Number of Plots   =   " << nplots << XMPI::endl;
        ss << "  List of Plots (PAR LAT LON REF SAMPLE R_MIN R_MAX):" << XMPI::endl;
        for (int i = 0; i < nplots; i++) ss << "    * " << splots[i]->verbose() << XMPI::endl;
        ss << "======================== Slice Plots =======================\n" << XMPI::endl;
        XMPI::cout << ss.str();
    }
}


void SlicePlot::plotUnweighted() const {
    
}

void SlicePlot::plotWeighted() const {
    
}


