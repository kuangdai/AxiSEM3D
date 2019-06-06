// PointwiseIONetCDF.cpp
// created by Kuangdai on 1-Jun-2017 
// NetCDF IO for point-wise receivers

#include "PointwiseIONetCDF.h"
#include "Parameters.h"
#include "NetCDF_Writer.h"
#include "NetCDF_Reader.h"
#include "XMPI.h"
#include <sstream>
#include <cstdio>
#include "PointwiseRecorder.h"
#include <fstream>

void PointwiseIONetCDF::initialize(int totalRecordSteps, int bufferSize, 
    const std::string &components, const std::vector<PointwiseInfo> &receivers,
    double srcLat, double srcLon, double srcDep) {
    mReceivers = &receivers;
    // source location
    mSrcLat = srcLat;
    mSrcLon = srcLon;
    mSrcDep = srcDep;
    // number
    int numRec = receivers.size();
    int numStrainRec = 0;
    for (const auto &rec: *mReceivers) {
        if (rec.mDumpStrain) {
            numStrainRec++;
        }
    }
    int numCurlRec = 0;
    for (const auto &rec: *mReceivers) {
        if (rec.mDumpCurl) {
            numCurlRec++;
        }
    }
    std::vector<int> allNumRec;
    XMPI::gather(numRec, allNumRec, false);
    if (XMPI::root()) {
        mMinRankWithRec = -1;
        for (int i = 0; i < allNumRec.size(); i++) {
            if (allNumRec[i] > 0) {
                mMinRankWithRec = i;
                break;
            }
        }
    }
    XMPI::bcast(mMinRankWithRec);
    if (mMinRankWithRec == -1) {
        // no receiver at all
        return;
    }
    
    // station names
    mVarNamesDisp.resize(numRec);
    mVarNamesStrain.resize(numStrainRec);
    mStrainIndex.resize(numStrainRec);
    mVarNamesCurl.resize(numCurlRec);
    mCurlIndex.resize(numCurlRec);
    std::vector<double> mylats, mylons, mydeps;
    int istrain = 0;
    int icurl = 0;
    for (int irec = 0; irec < numRec; irec++) {
        mVarNamesDisp[irec] = receivers[irec].mNetwork + "." + receivers[irec].mName;
        mVarNamesDisp[irec] += "." + components;
        if (receivers[irec].mDumpStrain) {
            mVarNamesStrain[istrain] = receivers[irec].mNetwork + "." + receivers[irec].mName;
            mVarNamesStrain[istrain] += ".RTZ.strain";
            mStrainIndex[istrain] = irec;
            istrain++;
        }
        if (receivers[irec].mDumpCurl) {
            mVarNamesCurl[icurl] = receivers[irec].mNetwork + "." + receivers[irec].mName;
            mVarNamesCurl[icurl] += ".RTZ.curl";
            mCurlIndex[icurl] = irec;
            icurl++;
        }
        mylats.push_back(receivers[irec].mLat);
        mylons.push_back(receivers[irec].mLon);
        mydeps.push_back(receivers[irec].mDep);
    }
    
    // dims
    std::vector<size_t> dimsTime;
    std::vector<size_t> dimsSeis;
    std::vector<size_t> dimsStrain;
    std::vector<size_t> dimsCurl;
    
    dimsTime.push_back(totalRecordSteps);
    dimsSeis.push_back(totalRecordSteps);
    dimsSeis.push_back(3);
    dimsStrain.push_back(totalRecordSteps);
    dimsStrain.push_back(6);
    dimsCurl.push_back(totalRecordSteps);
    dimsCurl.push_back(3);

    #ifndef _USE_PARALLEL_NETCDF
        // file
        int nfile = numRec / mMaxNumRecPerFile;
        if (nfile * mMaxNumRecPerFile < numRec) {
            nfile++;
        }
        for (int ifile = 0; ifile < nfile; ifile++) {
            mNetCDFs.push_back(new NetCDF_Writer());
        }
        if (mNetCDFs.size() > 1 && mAssemble) {
            throw std::runtime_error("PointwiseIONetCDF::initialize || mNetCDFs.size() > 1 && mAssemble. "
                                    " || Not implemented. Use netcdf_no_assemble");
        }
        
        // assemble
        if (!mAssemble) {
            std::vector<std::string> myRecKeys;
            for (int irec = 0; irec < numRec; irec++) {
                myRecKeys.push_back((*mReceivers)[irec].mNetwork + "." + (*mReceivers)[irec].mName);
            }
            std::vector<std::vector<std::string>> allRecKeys;
            XMPI::gather(myRecKeys, allRecKeys, false);
            if (XMPI::root()) {
                std::fstream fout(Parameters::sOutputDirectory + "/stations/rank_station.txt", std::fstream::out);
                fout << "# MPI_RANK NETWORK.NAME\n";
                for (int rank = 0; rank < XMPI::nproc(); rank++) {
                    for (int irec = 0; irec < allRecKeys[rank].size(); irec++) {
                        int ifile = irec / mMaxNumRecPerFile;
                        fout << rank + XMPI::nproc() * ifile << " " << allRecKeys[rank][irec] << "\n";
                    }
                }
            }
        }
    
        // open file on all ranks with receivers
        if (numRec == 0) {
            return;
        }
        
        for (int ifile = 0; ifile < nfile; ifile++) {
            std::stringstream fname;
            fname << Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.nc.rank" 
                << XMPI::rank() + XMPI::nproc() * ifile;
            mNetCDFs[ifile]->open(fname.str(), true);
            mNetCDFs[ifile]->defModeOn();
            // define time
            mNetCDFs[ifile]->defineVariable<double>("time_points", dimsTime);
            // define seismograms
            int istrain = 0, icurl = 0;
            for (int irec = 0; irec < numRec; irec++) {
                if (irec / mMaxNumRecPerFile == ifile) {
                    mNetCDFs[ifile]->defineVariable<Real>(mVarNamesDisp[irec], dimsSeis);
                    mNetCDFs[ifile]->addAttribute(mVarNamesDisp[irec], "latitude", mylats[irec]);
                    mNetCDFs[ifile]->addAttribute(mVarNamesDisp[irec], "longitude", mylons[irec]);
                    mNetCDFs[ifile]->addAttribute(mVarNamesDisp[irec], "depth", mydeps[irec]);
                    if ((*mReceivers)[irec].mDumpStrain) {
                        mNetCDFs[ifile]->defineVariable<Real>(mVarNamesStrain[istrain++], dimsStrain);
                    }
                    if ((*mReceivers)[irec].mDumpCurl) {
                        mNetCDFs[ifile]->defineVariable<Real>(mVarNamesCurl[icurl++], dimsCurl);
                    }
                } else {
                    if ((*mReceivers)[irec].mDumpStrain) {
                        istrain++;
                    }
                    if ((*mReceivers)[irec].mDumpCurl) {
                        icurl++;
                    }
                }
            }
            mNetCDFs[ifile]->defModeOff();
            // fill time with err values
            mNetCDFs[ifile]->fillConstant("time_points", dimsTime, NC_ERR_VALUE);
            // fill seismograms with err values
            istrain = 0;
            icurl = 0;
            for (int irec = 0; irec < numRec; irec++) {
                if (irec / mMaxNumRecPerFile == ifile) {
                    mNetCDFs[ifile]->fillConstant(mVarNamesDisp[irec], dimsSeis, (Real)NC_ERR_VALUE);
                    if ((*mReceivers)[irec].mDumpStrain) {
                        mNetCDFs[ifile]->fillConstant(mVarNamesStrain[istrain++], dimsStrain, (Real)NC_ERR_VALUE);
                    }
                    if ((*mReceivers)[irec].mDumpCurl) {
                        mNetCDFs[ifile]->fillConstant(mVarNamesCurl[icurl++], dimsCurl, (Real)NC_ERR_VALUE);
                    }
                } else {
                    if ((*mReceivers)[irec].mDumpStrain) {
                        istrain++;
                    }
                    if ((*mReceivers)[irec].mDumpCurl) {
                        icurl++;
                    }
                }
            }
            // source location
            mNetCDFs[ifile]->addAttribute("", "source_latitude", mSrcLat);
            mNetCDFs[ifile]->addAttribute("", "source_longitude", mSrcLon);
            mNetCDFs[ifile]->addAttribute("", "source_depth", mSrcDep);
            mNetCDFs[ifile]->flush();
        }
    #else
        // gather all station names 
        std::vector<std::vector<std::string>> allNamesDisp, allNamesStrain, allNamesCurl; 
        std::vector<std::vector<int>> allIndexStrain;
        std::vector<std::vector<int>> allIndexCurl;
        std::vector<std::vector<double>> allLats;
        std::vector<std::vector<double>> allLons;
        std::vector<std::vector<double>> allDeps;
        XMPI::gather(mVarNamesDisp, allNamesDisp, true);
        XMPI::gather(mVarNamesStrain, allNamesStrain, true);
        XMPI::gather(mStrainIndex, allIndexStrain, MPI_INT, true);
        XMPI::gather(mVarNamesCurl, allNamesCurl, true);
        XMPI::gather(mCurlIndex, allIndexCurl, MPI_INT, true);
        XMPI::gather(mylats, allLats, MPI_DOUBLE, true);
        XMPI::gather(mylons, allLons, MPI_DOUBLE, true);
        XMPI::gather(mydeps, allDeps, MPI_DOUBLE, true);
        
        mNetCDFs = std::vector<NetCDF_Writer *>(1, new NetCDF_Writer());
        // open file on min rank and define all variables
        std::string fname = Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.nc";
        if (XMPI::rank() == mMinRankWithRec) {
            mNetCDFs[0]->open(fname, true);
            mNetCDFs[0]->defModeOn();
            // define time
            mNetCDFs[0]->defineVariable<double>("time_points", dimsTime);
            // define seismograms
            for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
                for (int irec = 0; irec < allNamesDisp[iproc].size(); irec++) {
                    mNetCDFs[0]->defineVariable<Real>(allNamesDisp[iproc][irec], dimsSeis);
                    mNetCDFs[0]->addAttribute(allNamesDisp[iproc][irec], "latitude", allLats[iproc][irec]);
                    mNetCDFs[0]->addAttribute(allNamesDisp[iproc][irec], "longitude", allLons[iproc][irec]);
                    mNetCDFs[0]->addAttribute(allNamesDisp[iproc][irec], "depth", allDeps[iproc][irec]);
                }
                for (int irec = 0; irec < allNamesStrain[iproc].size(); irec++) {
                    mNetCDFs[0]->defineVariable<Real>(allNamesStrain[iproc][irec], dimsStrain);
                    mNetCDFs[0]->addAttribute(allNamesStrain[iproc][irec], "latitude", allLats[iproc][allIndexStrain[iproc][irec]]);
                    mNetCDFs[0]->addAttribute(allNamesStrain[iproc][irec], "longitude", allLons[iproc][allIndexStrain[iproc][irec]]);
                    mNetCDFs[0]->addAttribute(allNamesStrain[iproc][irec], "depth", allDeps[iproc][allIndexStrain[iproc][irec]]);
                }
                for (int irec = 0; irec < allNamesCurl[iproc].size(); irec++) {
                    mNetCDFs[0]->defineVariable<Real>(allNamesCurl[iproc][irec], dimsCurl);
                    mNetCDFs[0]->addAttribute(allNamesCurl[iproc][irec], "latitude", allLats[iproc][allIndexCurl[iproc][irec]]);
                    mNetCDFs[0]->addAttribute(allNamesCurl[iproc][irec], "longitude", allLons[iproc][allIndexCurl[iproc][irec]]);
                    mNetCDFs[0]->addAttribute(allNamesCurl[iproc][irec], "depth", allDeps[iproc][allIndexCurl[iproc][irec]]);
                }
            }
            mNetCDFs[0]->defModeOff();
            // fill time with err values
            mNetCDFs[0]->fillConstant("time_points", dimsTime, NC_ERR_VALUE);
            // fill seismograms with err values
            for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
                for (int irec = 0; irec < allNamesDisp[iproc].size(); irec++) {
                    mNetCDFs[0]->fillConstant<Real>(allNamesDisp[iproc][irec], dimsSeis, (Real)NC_ERR_VALUE);
                }
                for (int irec = 0; irec < allNamesStrain[iproc].size(); irec++) {
                    mNetCDFs[0]->fillConstant<Real>(allNamesStrain[iproc][irec], dimsStrain, (Real)NC_ERR_VALUE);
                }
                for (int irec = 0; irec < allNamesCurl[iproc].size(); irec++) {
                    mNetCDFs[0]->fillConstant<Real>(allNamesCurl[iproc][irec], dimsCurl, (Real)NC_ERR_VALUE);
                }
            }
            // source location
            mNetCDFs[0]->addAttribute("", "source_latitude", mSrcLat);
            mNetCDFs[0]->addAttribute("", "source_longitude", mSrcLon);
            mNetCDFs[0]->addAttribute("", "source_depth", mSrcDep);
            mNetCDFs[0]->close();
        }
        XMPI::barrier();
        mNetCDFs[0]->openParallel(fname);
    #endif
    
    // record postion in nc file
    mCurrentRow = 0;
}

void PointwiseIONetCDF::finalize() {
    if (mMinRankWithRec == -1) {
        // no receiver at all
        return;
    }

    // dispose writer
    for (int ifile = 0; ifile < mNetCDFs.size(); ifile++) {
        mNetCDFs[ifile]->close();
        delete mNetCDFs[ifile];
    }
    
    #ifdef _USE_PARALLEL_NETCDF
        return;
    #endif
    
    if (!mAssemble) {
        return;
    }
    
    // file name
    std::string oneFile = Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.nc";
    std::stringstream fname;
    fname << Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.nc.rank" << XMPI::rank();
    std::string locFile = fname.str();
    
    // merge
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(true);
    #endif
    
    // dims
    std::vector<size_t> dimsTime;
    std::vector<size_t> dimsSeis;
    std::vector<size_t> dimsStrain;
    std::vector<size_t> dimsCurl;
    dimsTime.push_back(mCurrentRow);
    dimsSeis.push_back(mCurrentRow);
    dimsSeis.push_back(3);
    dimsStrain.push_back(mCurrentRow);
    dimsStrain.push_back(6);
    dimsCurl.push_back(mCurrentRow);
    dimsCurl.push_back(3);
    
    // create file 
    if (XMPI::rank() == mMinRankWithRec) {
        // read time
        NetCDF_Reader nr;
        nr.open(locFile);
        RDColX times;
        nr.read1D("time_points", times);
        nr.close();
        
        // create file and write time
        NetCDF_Writer nw;
        nw.open(oneFile, true);
        nw.defModeOn();
        nw.defineVariable<double>("time_points", dimsTime);
        nw.defModeOff();
        nw.writeVariableWhole("time_points", times);
        
        // source location
        nw.addAttribute("", "source_latitude", mSrcLat);
        nw.addAttribute("", "source_longitude", mSrcLon);
        nw.addAttribute("", "source_depth", mSrcDep);
        nw.close();
    }
    XMPI::barrier();
    
    // write seismograms
    int numRec = mVarNamesDisp.size();
    int numStrainRec = mVarNamesStrain.size();
    int numCurlRec = mVarNamesCurl.size();
    for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
        if (iproc == XMPI::rank() && numRec > 0) {
            // open
            NetCDF_Reader nr;
            nr.open(locFile);
            NetCDF_Writer nw;
            nw.open(oneFile, false);
            
            // create variable
            nw.defModeOn();
            for (int irec = 0; irec < numRec; irec++) {
                nw.defineVariable<Real>(mVarNamesDisp[irec], dimsSeis);
                nw.addAttribute(mVarNamesDisp[irec], "latitude", (*mReceivers)[irec].mLat);
                nw.addAttribute(mVarNamesDisp[irec], "longitude", (*mReceivers)[irec].mLon);
                nw.addAttribute(mVarNamesDisp[irec], "depth", (*mReceivers)[irec].mDep);
            }
            for (int irec = 0; irec < numStrainRec; irec++) {
                nw.defineVariable<Real>(mVarNamesStrain[irec], dimsStrain);
                nw.addAttribute(mVarNamesStrain[irec], "latitude", (*mReceivers)[mStrainIndex[irec]].mLat);
                nw.addAttribute(mVarNamesStrain[irec], "longitude", (*mReceivers)[mStrainIndex[irec]].mLon);
                nw.addAttribute(mVarNamesStrain[irec], "depth", (*mReceivers)[mStrainIndex[irec]].mDep);
            }
            for (int irec = 0; irec < numCurlRec; irec++) {
                nw.defineVariable<Real>(mVarNamesCurl[irec], dimsCurl);
                nw.addAttribute(mVarNamesCurl[irec], "latitude", (*mReceivers)[mCurlIndex[irec]].mLat);
                nw.addAttribute(mVarNamesCurl[irec], "longitude", (*mReceivers)[mCurlIndex[irec]].mLon);
                nw.addAttribute(mVarNamesCurl[irec], "depth", (*mReceivers)[mCurlIndex[irec]].mDep);
            }
            nw.defModeOff();
            
            // read and write seismograms
            for (int irec = 0; irec < numRec; irec++) {
                RMatXX_RM seis;
                nr.read2D(mVarNamesDisp[irec], seis);
                nw.writeVariableWhole(mVarNamesDisp[irec], seis);
            }
            for (int irec = 0; irec < numStrainRec; irec++) {
                RMatXX_RM seis;
                nr.read2D(mVarNamesStrain[irec], seis);
                nw.writeVariableWhole(mVarNamesStrain[irec], seis);
            }
            for (int irec = 0; irec < numCurlRec; irec++) {
                RMatXX_RM seis;
                nr.read2D(mVarNamesCurl[irec], seis);
                nw.writeVariableWhole(mVarNamesCurl[irec], seis);
            }
            
            // close
            nw.close();
            nr.close();
        }
        XMPI::barrier();
    } 
    
    // delete local files
    if (numRec > 0) {
        std::remove(locFile.c_str());
    }
    
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(false);
    #endif
}

void PointwiseIONetCDF::dumpToFile(const RMatXX_RM &bufferDisp, 
    const RMatXX_RM &bufferStrain,
    const RMatXX_RM &bufferCurl,
    const RDColX &bufferTime, int bufferLine) {
    if (bufferLine == 0) {
        return;
    }
    
    // write time
    std::vector<size_t> start;
    std::vector<size_t> countDisp;
    std::vector<size_t> countStrain;    
    std::vector<size_t> countCurl;    
    start.push_back(mCurrentRow);
    countDisp.push_back(bufferLine);
    countStrain.push_back(bufferLine);
    countCurl.push_back(bufferLine);
    
    // update record postion in nc file
    mCurrentRow += bufferLine;
    
    int numRec = mVarNamesDisp.size();
    int numStrainRec = mVarNamesStrain.size();
    int numCurlRec = mVarNamesCurl.size();
    #ifndef _USE_PARALLEL_NETCDF
        if (numRec == 0) {
            return;
        }
        for (int ifile = 0; ifile < mNetCDFs.size(); ifile++) {
            mNetCDFs[ifile]->writeVariableChunk("time_points", 
                bufferTime.topRows(bufferLine), start, countDisp);
        }  
    #else
        if (XMPI::rank() == mMinRankWithRec) {
            mNetCDFs[0]->writeVariableChunk("time_points", 
                bufferTime.topRows(bufferLine), start, countDisp);
        }
    #endif
    
    // write seismograms
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(true);
    #endif
    start.push_back(0);
    countDisp.push_back(3);
    countStrain.push_back(6);
    countCurl.push_back(3);
    int istrain = 0, icurl = 0;
    for (int irec = 0; irec < numRec; irec++) {
        int ifile = irec / mMaxNumRecPerFile;
        mNetCDFs[ifile]->writeVariableChunk(mVarNamesDisp[irec], 
            bufferDisp.block(0, irec * 3, bufferLine, 3).eval(), 
            start, countDisp);
        if ((*mReceivers)[irec].mDumpStrain) {
            mNetCDFs[ifile]->writeVariableChunk(mVarNamesStrain[istrain++], 
                bufferStrain.block(0, irec * 6, bufferLine, 6).eval(), 
                start, countStrain);
        }
        if ((*mReceivers)[irec].mDumpCurl) {
            mNetCDFs[ifile]->writeVariableChunk(mVarNamesCurl[icurl++], 
                bufferCurl.block(0, irec * 3, bufferLine, 3).eval(), 
                start, countCurl);
        }
    }
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(false);
    #endif
    
    // flush to disk
    for (int ifile = 0; ifile < mNetCDFs.size(); ifile++) {
        mNetCDFs[ifile]->flush();
    }
}

