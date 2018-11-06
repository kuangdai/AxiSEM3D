// PointwiseIOAscii.cpp
// created by Kuangdai on 1-Jun-2017 
// ascii IO for point-wise receivers

#include "PointwiseIOAscii.h"
#include "Parameters.h"
#include "PointwiseRecorder.h"

void PointwiseIOAscii::initialize(int totalRecordSteps, int bufferSize, 
    const std::string &components, const std::vector<PointwiseInfo> &receivers,
    double srcLat, double srcLon, double srcDep) {
    // number
    int numRec = receivers.size();
    mFileNamesDisp.clear();    
    mFilesDisp.clear();
    mBufferDisp = RDMatXX::Zero(bufferSize, 4);
    mFileNamesStrain.clear();    
    mFilesStrain.clear();
    mBufferStrain = RDMatXX::Zero(bufferSize, 7);
    mFileNamesCurl.clear();    
    mFilesCurl.clear();
    mBufferCurl = RDMatXX::Zero(bufferSize, 4);
    
    // files
    std::string outdir = Parameters::sOutputDirectory + "/stations/";
    for (int irec = 0; irec < numRec; irec++) {
        std::string fname_base = outdir + receivers[irec].mNetwork + "." + receivers[irec].mName;
        std::string fname_disp = fname_base + "." + components + ".ascii";
        std::fstream *fs_disp = new std::fstream(fname_disp, std::fstream::out);
        if (!(*fs_disp)) {
            throw std::runtime_error("PointwiseIOAscii::initialize || "
                "Error opening ascii output file: || " + fname_disp
                + " || Use NetCDF instead of ascii if there are too many stations.");
        }
        mFileNamesDisp.push_back(fname_disp);
        mFilesDisp.push_back(fs_disp);
        
        // strain
        if (receivers[irec].mDumpStrain) {
            std::string fname_strain = fname_base + "." + "RTZ" + ".strain.ascii";
            std::fstream *fs_strain = new std::fstream(fname_strain, std::fstream::out);
            if (!(*fs_strain)) {
                throw std::runtime_error("PointwiseIOAscii::initialize || "
                    "Error opening ascii output file: || " + fname_strain
                    + " || Use NetCDF instead of ascii if there are too many stations.");
            }
            mFileNamesStrain.push_back(fname_strain);
            mFilesStrain.push_back(fs_strain);
        }
        
        // curl
        if (receivers[irec].mDumpCurl) {
            std::string fname_curl = fname_base + "." + "RTZ" + ".curl.ascii";
            std::fstream *fs_curl = new std::fstream(fname_curl, std::fstream::out);
            if (!(*fs_curl)) {
                throw std::runtime_error("PointwiseIOAscii::initialize || "
                    "Error opening ascii output file: || " + fname_curl
                    + " || Use NetCDF instead of ascii if there are too many stations.");
            }
            mFileNamesCurl.push_back(fname_curl);
            mFilesCurl.push_back(fs_curl);
        }
    }
}

void PointwiseIOAscii::finalize() {
    int numRec = mFileNamesDisp.size();
    for (int irec = 0; irec < numRec; irec++) {
        mFilesDisp[irec]->close();
        delete mFilesDisp[irec];
    }
    int numStrainRec = mFileNamesStrain.size();
    for (int irec = 0; irec < numStrainRec; irec++) {
        mFilesStrain[irec]->close();
        delete mFilesStrain[irec];
    }
    int numCurlRec = mFileNamesCurl.size();
    for (int irec = 0; irec < numCurlRec; irec++) {
        mFilesCurl[irec]->close();
        delete mFilesCurl[irec];
    }
}

void PointwiseIOAscii::dumpToFile(const RMatXX_RM &bufferDisp, 
    const RMatXX_RM &bufferStrain,
    const RMatXX_RM &bufferCurl,
    const RDColX &bufferTime, int bufferLine) {
    if (bufferLine == 0) {
        return;
    }
    
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(true);
    #endif
    
    int numRec = mFileNamesDisp.size();
    for (int irec = 0; irec < numRec; irec++) {
        mBufferDisp.topRows(bufferLine) << bufferTime.topRows(bufferLine), 
                                       bufferDisp.block(0, irec * 3, bufferLine, 3).cast<double>();
        (*mFilesDisp[irec]) << mBufferDisp.topRows(bufferLine).format(EIGEN_FMT) << std::endl;   
        mFilesDisp[irec]->flush();
    }
    
    int numStrainRec = mFileNamesStrain.size();
    for (int irec = 0; irec < numStrainRec; irec++) {
        mBufferStrain.topRows(bufferLine) << bufferTime.topRows(bufferLine), 
                                       bufferStrain.block(0, irec * 6, bufferLine, 6).cast<double>();
        (*mFilesStrain[irec]) << mBufferStrain.topRows(bufferLine).format(EIGEN_FMT) << std::endl;   
        mFilesStrain[irec]->flush();
    }
    
    int numCurlRec = mFileNamesCurl.size();
    for (int irec = 0; irec < numCurlRec; irec++) {
        mBufferCurl.topRows(bufferLine) << bufferTime.topRows(bufferLine), 
                                       bufferCurl.block(0, irec * 3, bufferLine, 3).cast<double>();
        (*mFilesCurl[irec]) << mBufferCurl.topRows(bufferLine).format(EIGEN_FMT) << std::endl;   
        mFilesCurl[irec]->flush();
    }
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(false);
    #endif
}

