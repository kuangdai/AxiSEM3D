// SlicePlot.cpp
// created by Kuangdai on 26-Nov-2016 
// plot a slice

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <sstream>
#include <fstream>

#include "SlicePlot.h"

#include "Geodesy.h"
#include "XMPI.h"
#include "Mesh.h"
#include "Parameters.h"
#include "ExodusModel.h"
#include "Quad.h"
#include "Domain.h"
#include "Element.h"

SlicePlot::SlicePlot(const std::string &params, const Mesh *mesh):
mMesh(mesh) {
    std::string mstr = params;
    std::vector<std::string> strs;
    boost::trim_if(mstr, boost::is_any_of("\t "));
    boost::split(strs, mstr, boost::is_any_of("$"), boost::token_compress_on);
    mParName = strs[0];
    int narg = numArgs(mParName);
    if (strs.size() < narg) throw std::runtime_error("SlicePlot::SlicePlot || "
        "Not enough arguments following ParSeries " + mParName + ".");
    try {
        if (narg >= 2) {
            // sample type
            if (boost::iequals(strs[1], "center")) {
                mSampleType = SampleTypes::Center;
            } else if (boost::iequals(strs[1], "vertex")) {
                mSampleType = SampleTypes::Vertex;
            } else if (boost::iequals(strs[1], "gllpnt")) {
                mSampleType = SampleTypes::GLLPnt;
            } else throw std::runtime_error("");
        }
        
        if (narg >= 4) {
            // phi
            mLat = boost::lexical_cast<double>(strs[2]);
            mLon = boost::lexical_cast<double>(strs[3]);
            RDCol3 rtpG;
            rtpG(0) = 1.;
            rtpG(1) = Geodesy::lat2Theta_d(mLat, 0.);
            rtpG(2) = Geodesy::lon2Phi(mLon);
            const RDCol3 &rtpS = Geodesy::rotateGlob2Src(rtpG, mMesh->mSrcLat, mMesh->mSrcLon, mMesh->mSrcDep);
            mPhi = rtpS(2);
        }
        
        if (narg >= 5) {
            // property type
            if (boost::iequals(strs[4], "1D")) {
                mRefType = PropertyRefTypes::Property1D;
            } else if (boost::iequals(strs[4], "3D")) {
                mRefType = PropertyRefTypes::Property3D;
            } else if (boost::iequals(strs[4], "perturb")) {
                mRefType = PropertyRefTypes::PropertyPerturb;
            } else throw std::runtime_error("");
        }
        
    } catch (std::exception) {
        throw std::runtime_error("SlicePlot::SlicePlot || "
            "Invalid parameter following ParSeries " + mParName + ".");
    }      
}

void SlicePlot::buildInparam(std::vector<SlicePlot *> &splots, 
    const Parameters &par, const Mesh *mesh, int verbose) {
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
    if (verbose && nplots > 0) {
        std::stringstream ss;
        ss << "\n======================== Slice Plots =======================" << XMPI::endl;
        ss << "  Number of Plots   =   " << nplots << XMPI::endl;
        ss << "  List of Plots (PAR SAMPLE LAT LON PHI REF):" << XMPI::endl;
        for (int i = 0; i < nplots; i++) ss << "    * " << splots[i]->verbose() << XMPI::endl;
        ss << "======================== Slice Plots =======================\n" << XMPI::endl;
        XMPI::cout << ss.str();
    }
}


void SlicePlot::plotUnweighted() const {
    plotPhysicalProperties();
    plotNu();
    plotRank(false);
}

void SlicePlot::plotWeighted() const {
    plotRank(true);
}

void SlicePlot::plotPhysicalProperties() const {
    if (!isPhysical(mParName)) return;
    
    // data size
    int ncol = 1;
    if (mSampleType == SampleTypes::Vertex) ncol = 4;
    if (mSampleType == SampleTypes::GLLPnt) ncol = nPntElem;
    int nrow = mMesh->mExModel->getNumQuads();
    RDMatXX data = RDMatXX::Zero(nrow, ncol);
    
    // gather local data
    for (int i = 0; i < mMesh->getNumQuads(); i++) {
        RDRowN quadData;
        if (boost::iequals(mParName, "undulation")) {
            quadData = mMesh->mQuads[i]->getUndulationOnSlice(mPhi);
        } else {
            quadData = mMesh->mQuads[i]->getMaterialOnSlice(mParName, mRefType, mPhi);
        }
        int row = mMesh->mQuads[i]->getQuadTag();
        if (mSampleType == SampleTypes::Center) {
            data(row, 0) = quadData(nPol / 2 * nPntEdge + nPol / 2);
        } else if (mSampleType == SampleTypes::Vertex) {
            data(row, 0) = quadData(0);
            data(row, 1) = quadData(nPol * nPntEdge);
            data(row, 2) = quadData(nPol * nPntEdge + nPol);
            data(row, 3) = quadData(nPol);
        } else {
            data.row(row) = quadData;
        }
    }
    
    // sum MPI
    XMPI::sumEigenDouble(data);
    
    // dump to file
    if (XMPI::root()) {
        std::string fname = Parameters::sOutputDirectory + "/plots/" + verbose('_') + ".txt";
        std::fstream fs(fname, std::fstream::out);
        fs << data << std::endl;
        fs.close();
    }
}

void SlicePlot::plotNu() const {
    if (!boost::iequals(mParName, "nu")) return;
    
    // data size
    int ncol = 1;
    if (mSampleType == SampleTypes::Vertex) ncol = 4;
    if (mSampleType == SampleTypes::GLLPnt) ncol = nPntElem;
    int nrow = mMesh->mExModel->getNumQuads();
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> data(nrow, ncol);
    data.setZero();
    
    // gather local data
    for (int i = 0; i < mMesh->getNumQuads(); i++) {
        int row = mMesh->mQuads[i]->getQuadTag();
        if (mSampleType == SampleTypes::Center) {
            data(row, 0) = mMesh->mQuads[i]->getNu(); // use max
        } else if (mSampleType == SampleTypes::Vertex) {
            data(row, 0) = mMesh->mQuads[i]->getPointNr(0, 0);
            data(row, 1) = mMesh->mQuads[i]->getPointNr(nPol, 0);
            data(row, 2) = mMesh->mQuads[i]->getPointNr(nPol, nPol);
            data(row, 3) = mMesh->mQuads[i]->getPointNr(0, nPol);
        } else {
            for (int ipol = 0; ipol <= nPol; ipol++) {
                for (int jpol = 0; jpol <= nPol; jpol++) {
                    int ipnt = ipol * nPntEdge + jpol;
                    data(row, ipnt) = mMesh->mQuads[i]->getPointNr(ipol, jpol);
                }
            }
        }
    }
    
    // sum MPI
    XMPI::sumEigenInt(data);
    
    // dump to file
    if (XMPI::root()) {
        std::string fname = Parameters::sOutputDirectory + "/plots/" + verbose('_') + ".txt";
        std::fstream fs(fname, std::fstream::out);
        fs << data << std::endl;
        fs.close();
    }
}

void SlicePlot::plotRank(bool weighted) const {
    bool plot = boost::iequals(mParName, "weighted") && weighted;
    plot = plot || (boost::iequals(mParName, "unweighted") && (!weighted));
    if (!plot) return;
    
    // data size
    int nrow = mMesh->mExModel->getNumQuads();
    IColX data = IColX::Zero(nrow);
    
    // gather local data
    for (int i = 0; i < mMesh->getNumQuads(); i++) {
        int row = mMesh->mQuads[i]->getQuadTag();
        data(row) = XMPI::rank();
    }
    
    // sum MPI
    XMPI::sumEigenInt(data);
    
    // dump to file
    if (XMPI::root()) {
        std::string fname = Parameters::sOutputDirectory + "/plots/" + verbose('_') + ".txt";
        std::fstream fs(fname, std::fstream::out);
        fs << data << std::endl;
        fs.close();
    }
}

void SlicePlot::plotEleType(const Domain &domain) const {
    if (!boost::iequals(mParName, "eleType")) return;
    
    // data size
    int nrow = mMesh->mExModel->getNumQuads();
    std::vector<std::string> data(nrow, "");
    std::vector<std::string> etype;
    std::vector<int> etag;
    
    // gather local data
    for (int iloc = 0; iloc < mMesh->getNumQuads(); iloc++) {
        int elemTag = mMesh->mQuads[iloc]->getElementTag();
        int quadTag = mMesh->mQuads[iloc]->getQuadTag();
        Element *elem = domain.getElement(elemTag);
        etype.push_back(elem->verbose());
        etag.push_back(quadTag);
    }
    
    // sum MPI
    std::vector<std::vector<std::string>> all_etype;
    XMPI::gather(etype, all_etype, false);
    std::vector<std::vector<int>> all_etag;
    XMPI::gather(etag, all_etag, MPI_INT, false);
    
    // dump to file
    if (XMPI::root()) {
        for (int i = 0; i < all_etype.size(); i++) {
            for (int j = 0; j < all_etype[i].size(); j++) {
                data[all_etag[i][j]] = all_etype[i][j];
            }
        }
        std::string fname = Parameters::sOutputDirectory + "/plots/" + verbose('_') + ".txt";
        std::fstream fs(fname, std::fstream::out);
        for (int i = 0; i < nrow; i++) fs << data[i] << std::endl;
        fs.close();
    }
}

void SlicePlot::plotMeasured(const RDColX &cost) const {
    if (!boost::iequals(mParName, "eleCost")) return;
    // dump to file
    if (XMPI::root()) {
        std::string fname = Parameters::sOutputDirectory + "/plots/" + verbose('_') + ".txt";
        std::fstream fs(fname, std::fstream::out);
        fs << cost << std::endl;
        fs.close();
    }
}

int SlicePlot::numArgs(const std::string &parName) const {
    if (boost::iequals(parName, "vp")) return 5;
    if (boost::iequals(parName, "vpv")) return 5;
    if (boost::iequals(parName, "vph")) return 5;
    if (boost::iequals(parName, "vs")) return 5;
    if (boost::iequals(parName, "vsv")) return 5;
    if (boost::iequals(parName, "vsh")) return 5;
    if (boost::iequals(parName, "rho")) return 5;
    if (boost::iequals(parName, "eta")) return 5;
    if (boost::iequals(parName, "Qmu")) return 5;
    if (boost::iequals(parName, "Qkappa")) return 5;
    if (boost::iequals(parName, "undulation")) return 4;
    if (boost::iequals(parName, "nu")) return 2;
    if (boost::iequals(parName, "eleType")) return 1;
    if (boost::iequals(parName, "eleCost")) return 1;
    if (boost::iequals(parName, "unweighted")) return 1;
    if (boost::iequals(parName, "weighted")) return 1;
    throw std::runtime_error("SlicePlot::numArgs || Unknown parameter name: " + parName + ".");
}

bool SlicePlot::isPhysical(const std::string &parName) const {
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
    throw std::runtime_error("SlicePlot::onSlice || Unknown parameter name: " + parName + ".");
}

std::string SlicePlot::verbose(char sp) const {
    std::stringstream ss;
    ss << mParName;
    int narg = numArgs(mParName);
    if (narg >= 2) {
        if (mSampleType == SampleTypes::Center) ss << sp << "center";
        if (mSampleType == SampleTypes::Vertex) ss << sp << "vertex";
        if (mSampleType == SampleTypes::GLLPnt) ss << sp << "gllpnt";
    } 
    if (narg >= 4) ss << sp << mLat << sp << mLon << sp << mPhi / degree;
    if (narg >= 5) {
        if (mRefType == PropertyRefTypes::Property1D) ss << sp << "1D";
        if (mRefType == PropertyRefTypes::Property3D) ss << sp << "3D";
        if (mRefType == PropertyRefTypes::PropertyPerturb) ss << sp << "perturb";
    }
    return ss.str();
}


