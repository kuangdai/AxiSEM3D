// Source.cpp
// created by Kuangdai on 8-May-2016
// base class of source
// we only consider point sources located on the axis

#include "Source.h"
#include "Quad.h"
#include "Domain.h"
#include "Element.h"
#include "SourceTerm.h"
#include "Mesh.h"
#include "XMath.h"
#include "SpectralConstants.h"
#include "XMPI.h"
#include "MultilevelTimer.h"

Source::Source(double depth, double lat, double lon):
mDepth(depth), mLatitude(lat), mLongitude(lon) {
    // handle singularity at poles
    if (std::abs(mLatitude - 90.) < tinyDouble) {
        mLatitude = 90.;
        mLongitude = 0.;
    }
    if (std::abs(mLatitude + 90.) < tinyDouble) {
        mLatitude = -90.;
        mLongitude = 0.;
    }
}

void Source::release(Domain &domain, const Mesh &mesh) const {
    MultilevelTimer::begin("Locate Source", 2);
    // locate local
    int myrank = XMPI::nproc();
    int locTag;
    RDColP interpFactZ;
    if (locate(mesh, locTag, interpFactZ)) {
        myrank = XMPI::rank();
    }

    // min recRank
    int myrank_min = XMPI::min(myrank);
    if (myrank_min == XMPI::nproc()) {
        throw std::runtime_error("Source::release || Error locating source.");
    }
    MultilevelTimer::end("Locate Source", 2);

    MultilevelTimer::begin("Compute Source", 2);
    // release to me
    if (myrank_min == XMPI::rank()) {
        // compute source term
        arPP_CMatX3 fouriers;
        const Quad *myQuad = mesh.getQuad(locTag);
        computeSourceFourier(*myQuad, interpFactZ, fouriers);
        // add to domain
        Element *myElem = domain.getElement(myQuad->getElementTag());
        domain.addSourceTerm(new SourceTerm(myElem, fouriers));
    }
    MultilevelTimer::end("Compute Source", 2);
}

bool Source::locate(const Mesh &mesh, int &locTag, RDColP &interpFactZ) const {
    MultilevelTimer::begin("R Source", 3);
    RDCol2 srcCrds = RDCol2::Zero();
    srcCrds(1) = mesh.computeRadiusRef(mDepth, mLatitude, mLongitude);
    MultilevelTimer::end("R Source", 3);

    // check range of subdomain
    if (srcCrds(0) > mesh.sMax() + tinySingle || srcCrds(0) < mesh.sMin() - tinySingle) {
        return false;
    }
    if (srcCrds(1) > mesh.zMax() + tinySingle || srcCrds(1) < mesh.zMin() - tinySingle) {
        return false;
    }
    // find host element
    RDCol2 srcXiEta;
    for (int iloc = 0; iloc < mesh.getNumQuads(); iloc++) {
        const Quad *quad = mesh.getQuad(iloc);
        if (!quad->isAxial() || quad->isFluid() || !quad->nearMe(srcCrds(0), srcCrds(1))) {
            continue;
        }
        if (quad->invMapping(srcCrds, srcXiEta)) {
            if (std::abs(srcXiEta(1)) <= 1.000001) {
                if (std::abs(srcXiEta(0) + 1.) > tinySingle) {
                    throw std::runtime_error("Source::locate || Bad source location.");
                }
                locTag = iloc;
                XMath::interpLagrange(srcXiEta(1), nPntEdge,
                    SpectralConstants::getP_GLL().data(), interpFactZ.data());
                return true;
            }
        }
    }
    return false;
}

#include "Parameters.h"
#include "Earthquake.h"
#include "PointForce.h"
#include "NullSource.h"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <cfloat>

void Source::parseLine(const std::string &line, const std::string &key, double &res) {
    std::stringstream ss;
    std::string keylong;
    ss << line;
    ss >> keylong;
    if (boost::icontains(keylong, key)) {
        ss >> res;
    }
}

void Source::checkValue(const std::string &key, double res) {
    if (res > DBL_MAX * .99) {
        throw std::runtime_error("Source::checkValue || "
            "Error initializing source parameter: " + key);
    }
}

void Source::buildInparam(Source *&src, const Parameters &par, int verbose) {
    if (src) {
        delete src;
    }

    // null source
    if (par.getValue<bool>("DEVELOP_NON_SOURCE_MODE")) {
        src = new NullSource();
        if (verbose) {
            XMPI::cout << src->verbose();
        }
        return;
    }

    std::string src_type = par.getValue<std::string>("SOURCE_TYPE");
    std::string src_file = par.getValue<std::string>("SOURCE_FILE");

    if (boost::iequals(src_type, "earthquake")) {
        std::string cmtfile = Parameters::sInputDirectory + "/" + src_file;
        double depth = DBL_MAX, lat = DBL_MAX, lon = DBL_MAX;
        double Mrr = DBL_MAX, Mtt = DBL_MAX, Mpp = DBL_MAX; 
        double Mrt = DBL_MAX, Mrp = DBL_MAX, Mtp = DBL_MAX;
        if (XMPI::root()) {
            std::fstream fs(cmtfile, std::fstream::in);
            if (!fs) {
                throw std::runtime_error("Source::buildInparam || "
                    "Error opening CMT data file: ||" + cmtfile);
            }
            std::string line;
            while (std::getline(fs, line)) {
                parseLine(line, "latitude", lat);
                parseLine(line, "longitude", lon);
                parseLine(line, "depth", depth);
                parseLine(line, "Mrr", Mrr);
                parseLine(line, "Mtt", Mtt);
                parseLine(line, "Mpp", Mpp);
                parseLine(line, "Mrt", Mrt);
                parseLine(line, "Mrp", Mrp);
                parseLine(line, "Mtp", Mtp);
            }
            checkValue("latitude", lat);
            checkValue("longitude", lon);
            checkValue("depth", depth);
            checkValue("Mrr", Mrr);
            checkValue("Mtt", Mtt);
            checkValue("Mpp", Mpp);
            checkValue("Mrt", Mrt);
            checkValue("Mrp", Mrp);
            checkValue("Mtp", Mtp);
            // unit
            depth *= 1e3;
            Mrr *= 1e-7;
            Mtt *= 1e-7;
            Mpp *= 1e-7;
            Mrt *= 1e-7;
            Mrp *= 1e-7;
            Mtp *= 1e-7;
            fs.close();
        }
        XMPI::bcast(depth);
        XMPI::bcast(lat);
        XMPI::bcast(lon);
        XMPI::bcast(Mrr);
        XMPI::bcast(Mtt);
        XMPI::bcast(Mpp);
        XMPI::bcast(Mrt);
        XMPI::bcast(Mrp);
        XMPI::bcast(Mtp);
        src = new Earthquake(depth, lat, lon, Mrr, Mtt, Mpp, Mrt, Mrp, Mtp);
    } else if (boost::iequals(src_type, "point_force")) {
        // point force
        std::string pointffile = Parameters::sInputDirectory + "/" + src_file;
        double depth = DBL_MAX, lat = DBL_MAX, lon = DBL_MAX;
        double f1 = DBL_MAX, f2 = DBL_MAX, f3 = DBL_MAX;
        if (XMPI::root()) {
            std::fstream fs(pointffile, std::fstream::in);
            if (!fs) {
                throw std::runtime_error("Source::buildInparam || "
                    "Error opening point force data file: ||" + pointffile);
            }
            std::string line;
            while (std::getline(fs, line)) {
                parseLine(line, "latitude", lat);
                parseLine(line, "longitude", lon);
                parseLine(line, "depth", depth);
                parseLine(line, "Ft", f1);
                parseLine(line, "Fp", f2);
                parseLine(line, "Fr", f3);
            }
            checkValue("latitude", lat);
            checkValue("longitude", lon);
            checkValue("depth", depth);
            checkValue("Ft", f1);
            checkValue("Fp", f2);
            checkValue("Fr", f3);
            // unit
            depth *= 1e3;
            fs.close();
        }
        XMPI::bcast(depth);
        XMPI::bcast(lat);
        XMPI::bcast(lon);
        XMPI::bcast(f1);
        XMPI::bcast(f2);
        XMPI::bcast(f3);
        src = new PointForce(depth, lat, lon, f1, f2, f3);
    } else {
        throw std::runtime_error("Source::buildInparam || Unknown source type: " + src_type);
    }
    
    if (verbose) {
        XMPI::cout << src->verbose();
    }
}
