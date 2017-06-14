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

void PointwiseIONetCDF::initialize(int totalRecordSteps, int bufferSize, bool ENZ,
    const std::vector<std::string> &names,
    const std::vector<std::string> &networks) {
    // number
    int numRec = names.size();
    if (numRec == 0) {
        return;
    }    
        
    // file
    mNetCDF = new NetCDF_Writer();
    
    // open file on all ranks
    std::stringstream fname;
    fname << Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.nc.rank" << XMPI::rank();
    mNetCDF->open(fname.str(), true);
    
    // define time variable
    std::vector<size_t> dims;
    dims.push_back(totalRecordSteps);
    mNetCDF->defineVariable("time_points", dims, (Real)-1.2345);
    
    // define seismogram variables
    dims.push_back(3);
    mVarNames.resize(numRec);
    for (int irec = 0; irec < numRec; irec++) {
        mVarNames[irec] = networks[irec] + "." + names[irec];
        mVarNames[irec] += ENZ ? ".ENZ" : ".RTZ";
        mNetCDF->defineVariable(mVarNames[irec], dims, (Real)-1.2345);
    }
    
    // record postion in nc file
    mCurrentRow = 0;
}

void PointwiseIONetCDF::finalize() {
    // dispose writer
    int numRec = mVarNames.size();
    if (numRec > 0) {
        mNetCDF->close();
        delete mNetCDF;
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
    
    int fileDefined = 0;
    for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
        if (iproc == XMPI::rank() && numRec > 0) {
            // reader
            NetCDF_Reader nr;
            nr.open(locFile);
            
            // init file
            if (!fileDefined) {
                // read time
                RColX times;
                nr.read1D("time_points", times);
                // create file
                NetCDF_Writer nw;
                nw.open(oneFile, true);
                // write time
                std::vector<size_t> dims;
                dims.push_back(times.rows());
                nw.defineVariable("time_points", dims, (Real)-1.2345);
                nw.writeVariableWhole("time_points", times);
                nw.close();
                // file defined
                fileDefined = true;
            }
            
            // read and write seismograms
            NetCDF_Writer nw;
            nw.open(oneFile, false);
            for (int irec = 0; irec < numRec; irec++) {
                // read seis
                RMatXX_RM seis;
                nr.read2D(mVarNames[irec], seis);
                // write seis
                std::vector<size_t> dims;
                dims.push_back(seis.rows());
                dims.push_back(seis.cols());
                nw.defineVariable(mVarNames[irec], dims, (Real)-1.2345);
                nw.writeVariableWhole(mVarNames[irec], seis);
            }
            nw.close();
            
            // close reader
            nr.close();
        }
        XMPI::bcast(fileDefined);
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
    mNetCDF->writeVariableChunk("time_points", 
        bufferTime.topRows(bufferLine), start, count);
    
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
    
    // record postion in nc file
    mCurrentRow += bufferLine;
}

