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
    std::vector<double> mylats, mylons, mydeps;
    int istrain = 0;
    for (int irec = 0; irec < numRec; irec++) {
        mVarNamesDisp[irec] = receivers[irec].mNetwork + "." + receivers[irec].mName;
        mVarNamesDisp[irec] += "." + components;
        if (receivers[irec].mDumpStrain) {
            mVarNamesStrain[istrain] = receivers[irec].mNetwork + "." + receivers[irec].mName;
            mVarNamesStrain[istrain] += ".RTZ.strain";
            mStrainIndex[istrain] = irec;
            istrain++;
        }
        mylats.push_back(receivers[irec].mLat);
        mylons.push_back(receivers[irec].mLon);
        mydeps.push_back(receivers[irec].mDep);
    }
    
    // dims
    std::vector<size_t> dimsTime;
    std::vector<size_t> dimsSeis;
    std::vector<size_t> dimsStrain;
    dimsTime.push_back(totalRecordSteps);
    dimsSeis.push_back(totalRecordSteps);
    dimsSeis.push_back(3);
    dimsStrain.push_back(totalRecordSteps);
    dimsStrain.push_back(6);

    // file
    mNetCDF = new NetCDF_Writer();
    #ifndef _USE_PARALLEL_NETCDF
        // open file on all ranks with receivers
        if (numRec == 0) {
            return;
        }
        std::stringstream fname;
        fname << Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.nc.rank" << XMPI::rank();
        mNetCDF->open(fname.str(), true);
        mNetCDF->defModeOn();
        // define time
        mNetCDF->defineVariable<Real>("time_points", dimsTime);
        // define seismograms
        for (int irec = 0; irec < numRec; irec++) {
            mNetCDF->defineVariable<Real>(mVarNamesDisp[irec], dimsSeis);
            mNetCDF->addAttribute(mVarNamesDisp[irec], "latitude", mylats[irec]);
            mNetCDF->addAttribute(mVarNamesDisp[irec], "longitude", mylons[irec]);
            mNetCDF->addAttribute(mVarNamesDisp[irec], "depth", mydeps[irec]);
        }
        for (int irec = 0; irec < numStrainRec; irec++) {
            mNetCDF->defineVariable<Real>(mVarNamesStrain[irec], dimsStrain);
            mNetCDF->addAttribute(mVarNamesStrain[irec], "latitude", mylats[mStrainIndex[irec]]);
            mNetCDF->addAttribute(mVarNamesStrain[irec], "longitude", mylons[mStrainIndex[irec]]);
            mNetCDF->addAttribute(mVarNamesStrain[irec], "depth", mydeps[mStrainIndex[irec]]);
        }
        mNetCDF->defModeOff();
        // fill time with err values
        mNetCDF->fillConstant("time_points", dimsTime, (Real)NC_ERR_VALUE);
        // fill seismograms with err values
        for (int irec = 0; irec < numRec; irec++) {
            mNetCDF->fillConstant(mVarNamesDisp[irec], dimsSeis, (Real)NC_ERR_VALUE);
        }
        for (int irec = 0; irec < numStrainRec; irec++) {
            mNetCDF->fillConstant(mVarNamesStrain[irec], dimsStrain, (Real)NC_ERR_VALUE);
        }
        // source location
        mNetCDF->addAttribute("", "source_latitude", mSrcLat);
        mNetCDF->addAttribute("", "source_longitude", mSrcLon);
        mNetCDF->addAttribute("", "source_depth", mSrcDep);
        mNetCDF->flush();
    #else
        // gather all station names 
        std::vector<std::vector<std::string>> allNamesDisp, allNamesStrain; 
        std::vector<std::vector<int>> allIndexStrain;
        std::vector<std::vector<double>> allLats;
        std::vector<std::vector<double>> allLons;
        std::vector<std::vector<double>> allDeps;
        XMPI::gather(mVarNamesDisp, allNamesDisp, true);
        XMPI::gather(mVarNamesStrain, allNamesStrain, true);
        XMPI::gather(mStrainIndex, allIndexStrain, MPI_INT, true);
        XMPI::gather(mylats, allLats, MPI_DOUBLE, true);
        XMPI::gather(mylons, allLons, MPI_DOUBLE, true);
        XMPI::gather(mydeps, allDeps, MPI_DOUBLE, true);
        
        // open file on min rank and define all variables
        std::string fname = Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.nc";
        if (XMPI::rank() == mMinRankWithRec) {
            mNetCDF->open(fname, true);
            mNetCDF->defModeOn();
            // define time
            mNetCDF->defineVariable<Real>("time_points", dimsTime);
            // define seismograms
            for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
                for (int irec = 0; irec < allNamesDisp[iproc].size(); irec++) {
                    mNetCDF->defineVariable<Real>(allNamesDisp[iproc][irec], dimsSeis);
                    mNetCDF->addAttribute(allNamesDisp[iproc][irec], "latitude", allLats[iproc][irec]);
                    mNetCDF->addAttribute(allNamesDisp[iproc][irec], "longitude", allLons[iproc][irec]);
                    mNetCDF->addAttribute(allNamesDisp[iproc][irec], "depth", allDeps[iproc][irec]);
                }
                for (int irec = 0; irec < allNamesStrain[iproc].size(); irec++) {
                    mNetCDF->defineVariable<Real>(allNamesStrain[iproc][irec], dimsStrain);
                    mNetCDF->addAttribute(allNamesStrain[iproc][irec], "latitude", allLats[iproc][allIndexStrain[iproc][irec]]);
                    mNetCDF->addAttribute(allNamesStrain[iproc][irec], "longitude", allLons[iproc][allIndexStrain[iproc][irec]]);
                    mNetCDF->addAttribute(allNamesStrain[iproc][irec], "depth", allDeps[iproc][allIndexStrain[iproc][irec]]);
                }
            }
            mNetCDF->defModeOff();
            // fill time with err values
            mNetCDF->fillConstant("time_points", dimsTime, (Real)NC_ERR_VALUE);
            // fill seismograms with err values
            for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
                for (int irec = 0; irec < allNamesDisp[iproc].size(); irec++) {
                    mNetCDF->fillConstant<Real>(allNamesDisp[iproc][irec], dimsSeis, (Real)NC_ERR_VALUE);
                }
                for (int irec = 0; irec < allNamesStrain[iproc].size(); irec++) {
                    mNetCDF->fillConstant<Real>(allNamesStrain[iproc][irec], dimsStrain, (Real)NC_ERR_VALUE);
                }
            }
            // source location
            mNetCDF->addAttribute("", "source_latitude", mSrcLat);
            mNetCDF->addAttribute("", "source_longitude", mSrcLon);
            mNetCDF->addAttribute("", "source_depth", mSrcDep);
            mNetCDF->close();
        }
        XMPI::barrier();
        mNetCDF->openParallel(fname);
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
    mNetCDF->close();
    delete mNetCDF;
    
    #ifdef _USE_PARALLEL_NETCDF
        return;
    #endif
    
    if (!mAssemble) {
        int numRec = mVarNamesDisp.size();
        std::vector<std::string> myRecKeys;
        for (int irec = 0; irec < numRec; irec++) {
            myRecKeys.push_back((*mReceivers)[irec].mNetwork + "." + (*mReceivers)[irec].mName);
        }
        std::vector<std::vector<std::string>> allRecKeys;
        XMPI::gather(myRecKeys, allRecKeys, false);
        if (XMPI::root()) {
            std::fstream fout(Parameters::sOutputDirectory + "/stations/station_rank.txt", std::fstream::out);
            for (int rank = 0; rank < XMPI::nproc(); rank++) {
                if (allRecKeys[rank].size() > 0) {
                    fout << "RANK " << rank << ":\n";
                    for (int irec = 0; irec < allRecKeys[rank].size(); irec++) {
                        fout << allRecKeys[rank][irec] << "\n";
                    }
                    fout << "\n";    
                }
            }
        }
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
    dimsTime.push_back(mCurrentRow);
    dimsSeis.push_back(mCurrentRow);
    dimsSeis.push_back(3);
    dimsStrain.push_back(mCurrentRow);
    dimsStrain.push_back(6);
    
    // create file 
    if (XMPI::rank() == mMinRankWithRec) {
        // read time
        NetCDF_Reader nr;
        nr.open(locFile);
        RColX times;
        nr.read1D("time_points", times);
        nr.close();
        
        // create file and write time
        NetCDF_Writer nw;
        nw.open(oneFile, true);
        nw.defModeOn();
        nw.defineVariable<Real>("time_points", dimsTime);
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

void PointwiseIONetCDF::dumpToFile(const RMatXX_RM &bufferDisp, const RMatXX_RM &bufferStrain,
    const RColX &bufferTime, int bufferLine) {
    if (bufferLine == 0) {
        return;
    }
    
    // write time
    std::vector<size_t> start;
    std::vector<size_t> countDisp;
    std::vector<size_t> countStrain;    
    start.push_back(mCurrentRow);
    countDisp.push_back(bufferLine);
    countStrain.push_back(bufferLine);
    
    // update record postion in nc file
    mCurrentRow += bufferLine;
    
    int numRec = mVarNamesDisp.size();
    int numStrainRec = mVarNamesStrain.size();
    #ifndef _USE_PARALLEL_NETCDF
        if (numRec == 0) {
            return;
        }  
        mNetCDF->writeVariableChunk("time_points", 
            bufferTime.topRows(bufferLine), start, countDisp);
    #else
        if (XMPI::rank() == mMinRankWithRec) {
            mNetCDF->writeVariableChunk("time_points", 
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
    for (int irec = 0; irec < numRec; irec++) {
        mNetCDF->writeVariableChunk(mVarNamesDisp[irec], 
            bufferDisp.block(0, irec * 3, bufferLine, 3).eval(), 
            start, countDisp);
    }
    for (int irec = 0; irec < numStrainRec; irec++) {
        mNetCDF->writeVariableChunk(mVarNamesStrain[irec], 
            bufferStrain.block(0, irec * 6, bufferLine, 6).eval(), 
            start, countStrain);
    }
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(false);
    #endif
    
    // flush to disk
    mNetCDF->flush();
}

