// PointwiseIOASDF.cpp
// created by Kuangdai on 1-Jun-2017 
// NetCDF IO for point-wise receivers

#include "PointwiseIOASDF.h"
#include "Parameters.h"
#include "NetCDF_Writer.h"
#include "NetCDF_Reader.h"
#include "XMPI.h"
#include <sstream>
#include <cstdio>
#include "PointwiseRecorder.h"

void PointwiseIOASDF::initialize(int totalRecordSteps, int bufferSize, bool ENZ,
    const std::vector<PointwiseInfo> &receivers) {
    // number
    int numRec = receivers.size();
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
        mVarNames[irec] = receivers[irec].mNetwork + "." + receivers[irec].mName;
        mVarNames[irec] += ENZ ? ".ENZ" : ".RTZ";
        mNetCDF->defineVariable(mVarNames[irec], dims, (Real)-1.2345);
        // get info for station XML
        mNetworks.push_back(receivers[irec].mNetwork);
        mNames.push_back(receivers[irec].mName);
        mLats.push_back(receivers[irec].mLat);
        mLons.push_back(receivers[irec].mLon);
        mDeps.push_back(receivers[irec].mDep);
    }
    mENZ = ENZ;
    
    // record postion in nc file
    mCurrentRow = 0;
}

void PointwiseIOASDF::finalize() {
    // dispose writer
    int numRec = mVarNames.size();
    if (numRec > 0) {
        mNetCDF->close();
        delete mNetCDF;
    }
    
    // file name
    std::string oneFile = Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.asdf";
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
                // create file
                NetCDF_Writer nw;
                nw.open(oneFile, true);
                
                // create groups
                nw.createGroup("AuxiliaryData");
                nw.createGroup("Provenance");
                nw.createGroup("Waveforms");
                
                // create QuakeML
                createQuakeML(nw);
                
                // close
                nw.close();
                // file defined
                fileDefined = true;
            }
            
            // read and write seismograms
            NetCDF_Writer nw;
            nw.open(oneFile, false);
            for (int irec = 0; irec < numRec; irec++) {
                nw.goToFileRoot();
                std::string key = mNetworks[irec] + "." + mNames[irec];
                
                // create group
                nw.createGroup(key);
                nw.goToGroup(key);
                
                // create station XML
                createStationML(nw, irec);
                
                // read seis
                RMatXX_RM seis;
                nr.read2D(mVarNames[irec], seis);
                // write seis
                std::vector<size_t> dims;
                dims.push_back(seis.rows());
                std::string comp = mENZ? "ENZ" : "RTZ";
                for (int i = 0; i < 3; i++) {
                    std::string varName = key + "." + comp.c_str()[i];
                    nw.defineVariable(varName, dims, (Real)-1.2345);
                    nw.writeVariableWhole(varName, seis.col(i).eval());
                }
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

void PointwiseIOASDF::dumpToFile(const RMatXX_RM &bufferDisp, 
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

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
void PointwiseIOASDF::createQuakeML(NetCDF_Writer &nw) {
    std::string quakeStr;
    std::fstream fs(Parameters::sInputDirectory + "/ASDF/quake.xml");
    if (fs) {
        // user provided
        quakeStr = std::string((std::istreambuf_iterator<char>(fs)), 
            std::istreambuf_iterator<char>());
    } else {
        // built-in
        std::string path = projectDirectory + "/src/core/output/pointwise/minimum_quake.xml";
        fs.open(path);
        if (!fs) {
            throw std::runtime_error("PointwiseIOASDF::createQuakeML || "
                "Error opening minimum_quake.xml at directory: ||" + path);
        }
        quakeStr = std::string((std::istreambuf_iterator<char>(fs)), 
            std::istreambuf_iterator<char>());
        boost::replace_first(quakeStr, "@LAT@", boost::lexical_cast<std::string>(mSrcLat));
        boost::replace_first(quakeStr, "@LON@", boost::lexical_cast<std::string>(mSrcLon));
        boost::replace_first(quakeStr, "@DEP@", boost::lexical_cast<std::string>(mSrcDep));
    }
    fs.close();
    nw.writeStringInByte("QuakeML", quakeStr);
}

void PointwiseIOASDF::createStationML(NetCDF_Writer &nw, int irec) {
    std::string stationStr;
    std::string key = mNetworks[irec] + "." + mNames[irec];
    std::fstream fs(Parameters::sInputDirectory + "/ASDF/" + key + ".xml");
    if (fs) {
        // user provided
        stationStr = std::string((std::istreambuf_iterator<char>(fs)), 
            std::istreambuf_iterator<char>());
    } else {
        // built-in
        std::string path = projectDirectory + "/src/core/output/pointwise/minimum_station.xml";
        fs.open(path);
        if (!fs) {
            throw std::runtime_error("PointwiseIOASDF::createStationML || "
                "Error opening minimum_station.xml at directory: ||" + path);
        }
        stationStr = std::string((std::istreambuf_iterator<char>(fs)), 
            std::istreambuf_iterator<char>());
        boost::replace_first(stationStr, "@NETWORK@", boost::lexical_cast<std::string>(mNetworks[irec]));
        boost::replace_first(stationStr, "@NAME@", boost::lexical_cast<std::string>(mNames[irec]));
        boost::replace_first(stationStr, "@LAT@", boost::lexical_cast<std::string>(mLats[irec]));
        boost::replace_first(stationStr, "@LON@", boost::lexical_cast<std::string>(mLons[irec]));
        boost::replace_first(stationStr, "@DEP@", boost::lexical_cast<std::string>(mDeps[irec]));
    }
    fs.close();
    nw.writeStringInByte("StationXML", stationStr);
}
