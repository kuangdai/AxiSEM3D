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
#include "Mapping.h"
#include "XMath.h"
#include "SpectralConstants.h"
#include "XMPI.h"
#include "pointf.h"
#include "XTimer.h"

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
    XTimer::begin("Locate Source", 2);
    // locate local
    int myrank = XMPI::nproc();
    int locTag;
    RDColP interpFactZ;
    if (locate(mesh, locTag, interpFactZ)) myrank = XMPI::rank();

    // min recRank
    int myrank_min = XMPI::min(myrank);
    if (myrank_min == XMPI::nproc()) {
        throw std::runtime_error("Source::release || Error locating source.");
    }
    XTimer::end("Locate Source", 2);

    XTimer::begin("Compute Source", 2);
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
    XTimer::end("Compute Source", 2);
}

bool Source::locate(const Mesh &mesh, int &locTag, RDColP &interpFactZ) const {
    XTimer::begin("R Source", 3);
    RDCol2 srcCrds = RDCol2::Zero();
    srcCrds(1) = mesh.computeRadiusRef(mDepth, mLatitude, mLongitude);
    XTimer::end("R Source", 3);

    if (srcCrds(0) > mesh.sMax() + tinySingle || srcCrds(0) < mesh.sMin() - tinySingle) return false;
    if (srcCrds(1) > mesh.zMax() + tinySingle || srcCrds(1) < mesh.zMin() - tinySingle) return false;
    RDCol2 srcXiEta;
    for (int iloc = 0; iloc < mesh.getNumQuads(); iloc++) {
        const Quad *quad = mesh.getQuad(iloc);
        if (!quad->isAxial() || quad->isFluid()) continue;
        if (!quad->nearMe(srcCrds(0), srcCrds(1))) continue;
        if (quad->invMapping(srcCrds, srcXiEta)) {
            if (std::abs(srcXiEta(1)) <= 1.000001) {
                if (std::abs(srcXiEta(0) + 1.) > tinySingle)
                    throw std::runtime_error("Source::locate || Bad source location.");
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
#include "NullSource.h"
#include <fstream>
#include <boost/algorithm/string.hpp>

void Source::buildInparam(Source *&src, const Parameters &par, int verbose) {
    if (src) delete src;

    // null source
    if (par.getValue<bool>("DEVELOP_NON_SOURCE_MODE")) {
        src = new NullSource();
        if (verbose) XMPI::cout << src->verbose();
        return;
    }

    std::string src_type = par.getValue<std::string>("SOURCE_TYPE");
    std::string src_file = par.getValue<std::string>("SOURCE_FILE");

    if (boost::iequals(src_type, "earthquake")) {
        std::string cmtfile = Parameters::sInputDirectory + "/" + src_file;
        double depth, lat, lon;
        double Mrr, Mtt, Mpp, Mrt, Mrp, Mtp;
        if (XMPI::root()) {
            std::fstream fs(cmtfile, std::fstream::in);
            if (!fs) throw std::runtime_error("Earthquake::Earthquake || "
                "Error opening CMT data file: ||" + cmtfile);
            std::string junk;
            std::getline(fs, junk);
            std::getline(fs, junk);
            std::getline(fs, junk);
            std::getline(fs, junk);
            fs >> junk >> lat;
            fs >> junk >> lon;
            fs >> junk >> depth;
            fs >> junk >> Mrr;
            fs >> junk >> Mtt;
            fs >> junk >> Mpp;
            fs >> junk >> Mrt;
            fs >> junk >> Mrp;
            fs >> junk >> Mtp;
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
        if (verbose) XMPI::cout << src->verbose();
    } else {
        // point force
        std::string pointffile = Parameters::sInputDirectory + "/" + src_file;
  double depth, lat, lon;
  double f1, f2, f3;
  if (XMPI::root()) {
      std::fstream fs(pointffile, std::fstream::in);
      if (!fs) throw std::runtime_error("pointf::pointf || "
          "Error opening point force data file: ||" + pointffile);
      std::string junk;
      std::getline(fs, junk);
      std::getline(fs, junk);
      std::getline(fs, junk);
      std::getline(fs, junk);
      fs >> junk >> lat;
      fs >> junk >> lon;
      fs >> junk >> depth;
      fs >> junk >> f1;
      fs >> junk >> f2;
      fs >> junk >> f3;
      depth *= 1e3;
      fs.close();
  }
  XMPI::bcast(depth);
  XMPI::bcast(lat);
  XMPI::bcast(lon);
  XMPI::bcast(f1);
  XMPI::bcast(f2);
  XMPI::bcast(f3);
  src = new pointf(depth, lat, lon, f1,f2,f3);
  if (verbose) XMPI::cout << src->verbose();

    }
}
