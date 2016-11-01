// STF.cpp
// created by Kuangdai on 9-May-2016 
// source time function

#include "STF.h"
#include "SourceTimeFunction.h"
#include "Domain.h"

void STF::release(Domain &domain) const {
    std::vector<Real> ts(mSTF.begin(), mSTF.end());
    domain.setSTF(new SourceTimeFunction(ts, mDeltaT, mShift)); 
}

#include "XMPI.h"
#include "Parameters.h"
#include "ErfSTF.h"
#include <fstream>

void STF::buildInparam(STF *&stf, const Parameters &par, double dt, int verbose) {
    if (stf) delete stf;
    // max total steps
    int maxTotalSteps = INT_MAX;
    int enforceMaxSteps = par.getValue<int>("DEVELOP_MAX_TIME_STEPS");
    if (enforceMaxSteps > 0) maxTotalSteps = enforceMaxSteps;
    
    // read half duration
    std::string cmtfile = Parameters::sInputDirectory + "/CMTSOLUTION";
    double hdur;
    if (XMPI::root()) {
        std::fstream fs(cmtfile, std::fstream::in);
        if (!fs) throw std::runtime_error("STF::buildInparam || "
            "Error opening CMTSOLUTION data file: ||" + cmtfile);
        std::string line;
        while (getline(fs, line)) {
            try {
                std::vector<std::string> strs;
                boost::trim_if(line, boost::is_any_of("\t "));
                boost::split(strs, line, boost::is_any_of("\t "), boost::token_compress_on);
                if (boost::iequals(strs[0], "half")) {
                    hdur = boost::lexical_cast<double>(strs[strs.size() - 1]);
                    break;
                }
            } catch(std::exception) {
                throw std::runtime_error("STF::buildInparam || "
                    "Error reading half duration in CMTSOLUTION data file: ||" + cmtfile);
            }
        }
        fs.close();
    }
    XMPI::bcast(hdur);
    // Heaviside
    if (hdur < 5. * dt) hdur = 5. * dt;
    double decay = 1.628;
    double duration = par.getValue<double>("TIME_RECORD_LENGTH");
    stf = new ErfSTF(dt, duration, hdur, decay, maxTotalSteps);
    // verbose 
    if (verbose == 2) XMPI::cout << stf->verbose();
}

