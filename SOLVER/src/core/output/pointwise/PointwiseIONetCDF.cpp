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

void PointwiseIONetCDF::initialize(int totalRecordSteps, int bufferSize, bool ENZ,
    const std::vector<PointwiseInfo> &receivers) {
    // number
    int numRec = receivers.size();
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
    mVarNames.resize(numRec);
    for (int irec = 0; irec < numRec; irec++) {
        mVarNames[irec] = receivers[irec].mNetwork + "." + receivers[irec].mName;
        mVarNames[irec] += ENZ ? ".ENZ" : ".RTZ";
    }
    
    // dims
    std::vector<size_t> dimsTime;
    std::vector<size_t> dimsSeis;
    dimsTime.push_back(totalRecordSteps);
    dimsSeis.push_back(totalRecordSteps);
    dimsSeis.push_back(3);

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
            mNetCDF->defineVariable<Real>(mVarNames[irec], dimsSeis);
        }
        mNetCDF->defModeOff();
        // fill time with err values
        mNetCDF->fillConstant("time_points", dimsTime, (Real)NC_ERR_VALUE);
        // fill seismograms with err values
        for (int irec = 0; irec < numRec; irec++) {
            mNetCDF->fillConstant(mVarNames[irec], dimsSeis, (Real)NC_ERR_VALUE);
        }
        mNetCDF->flush();
    #else
        // gather all station names 
        std::vector<std::vector<std::string>> allNames; 
        XMPI::gather(mVarNames, allNames, true);
        // open file on min rank and define all variables
        std::string fname = Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.nc";
        if (XMPI::rank() == mMinRankWithRec) {
            mNetCDF->open(fname, true);
            mNetCDF->defModeOn();
            // define time
            mNetCDF->defineVariable<Real>("time_points", dimsTime);
            // define seismograms
            for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
                for (int irec = 0; irec < allNames[iproc].size(); irec++) {
                    mNetCDF->defineVariable<Real>(allNames[iproc][irec], dimsSeis);
                }
            }
            mNetCDF->defModeOff();
            // fill time with err values
            mNetCDF->fillConstant("time_points", dimsTime, (Real)NC_ERR_VALUE);
            // fill seismograms with err values
            for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
                for (int irec = 0; irec < allNames[iproc].size(); irec++) {
                    mNetCDF->fillConstant<Real>(allNames[iproc][irec], dimsSeis, (Real)NC_ERR_VALUE);
                }
            }
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
    
    // file name
    int numRec = mVarNames.size();
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
    dimsTime.push_back(mCurrentRow);
    dimsSeis.push_back(mCurrentRow);
    dimsSeis.push_back(3);
    
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
        nw.close();
    }
    XMPI::barrier();
    
    // write seismograms
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
                nw.defineVariable<Real>(mVarNames[irec], dimsSeis);
            }
            nw.defModeOff();
            
            // read and write seismograms
            for (int irec = 0; irec < numRec; irec++) {
                RMatXX_RM seis;
                nr.read2D(mVarNames[irec], seis);
                nw.writeVariableWhole(mVarNames[irec], seis);
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
    const RColX &bufferTime, int bufferLine) {
    int numRec = mVarNames.size();
    if (numRec == 0) {
        return;
    }  
    if (bufferLine == 0) {
        return;
    }
    
    // write time
    std::vector<size_t> start;
    std::vector<size_t> count;
    start.push_back(mCurrentRow);
    count.push_back(bufferLine);
    #ifndef _USE_PARALLEL_NETCDF
        mNetCDF->writeVariableChunk("time_points", 
            bufferTime.topRows(bufferLine), start, count);
    #else
        if (XMPI::rank() == mMinRankWithRec) {
            mNetCDF->writeVariableChunk("time_points", 
                bufferTime.topRows(bufferLine), start, count);
        }
    #endif
    
    // write seismograms
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(true);
    #endif
    start.push_back(0);
    count.push_back(3);
    for (int irec = 0; irec < numRec; irec++) {
        mNetCDF->writeVariableChunk(mVarNames[irec], 
            bufferDisp.block(0, irec * 3, bufferLine, 3).eval(), 
            start, count);
    }
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(false);
    #endif
    
    // flush to disk
    mNetCDF->flush();
    
    // record postion in nc file
    mCurrentRow += bufferLine;
}

