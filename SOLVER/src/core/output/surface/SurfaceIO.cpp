// SurfaceIO.cpp
// created by Kuangdai on 28-Nov-2017 
// NetCDF IO for surface wavefield

#include "SurfaceIO.h"
#include "Parameters.h"
#include "NetCDF_Writer.h"
#include "NetCDF_Reader.h"
#include "XMPI.h"
#include <sstream>
#include <cstdio>
#include "SurfaceInfo.h"
#include "eigenp.h"
#include "Geodesy.h"
#include "SpectralConstants.h"
#include <fstream>

void SurfaceIO::initialize(int totalRecordSteps, int bufferSize,
    const std::vector<SurfaceInfo> &surfaceInfo,
    double srcLat, double srcLon, double srcDep) {
    // source location
    mSrcLat = srcLat;
    mSrcLon = srcLon;
    mSrcDep = srcDep;
    // number
    int numEle = surfaceInfo.size();
    std::vector<int> allNumEle;
    XMPI::gather(numEle, allNumEle, false);
    if (XMPI::root()) {
        mMinRankWithEle = -1;
        for (int i = 0; i < allNumEle.size(); i++) {
            if (allNumEle[i] > 0) {
                mMinRankWithEle = i;
                break;
            }
        }
    }
    XMPI::bcast(mMinRankWithEle);
    if (mMinRankWithEle == -1) {
        // no element at all
        return;
    }
    
    // nc variable names
    mVarNames.resize(numEle);
    mNu.resize(numEle);
    for (int iele = 0; iele < numEle; iele++) {
        std::stringstream ss;
        ss << "edge_" << surfaceInfo[iele].getGlobalTag();
        mVarNames[iele] = ss.str();
        mNu[iele] = surfaceInfo[iele].getMaxNu();
    }
    
    // theta
    int numEleGlob = XMPI::sum(numEle);
    RDMatXX_RM theta = RDMatXX::Zero(numEleGlob, 2);
    for (int iele = 0; iele < numEle; iele++) {
        theta(surfaceInfo[iele].getGlobalTag(), 0) = surfaceInfo[iele].getTheta0();
        theta(surfaceInfo[iele].getGlobalTag(), 1) = surfaceInfo[iele].getTheta1();
    }
    XMPI::sumEigenDouble(theta);

    // dims
    std::vector<size_t> dimsTime;
    std::vector<size_t> dimsSeis;
    std::vector<size_t> dimsTheta;
    std::vector<size_t> dimsGLL;
    dimsTime.push_back(totalRecordSteps);
    dimsSeis.push_back(totalRecordSteps);
    dimsSeis.push_back(0);
    dimsTheta.push_back(numEleGlob);
    dimsTheta.push_back(2);
    dimsGLL.push_back(nPntEdge);
    
    // file
    mNetCDF = new NetCDF_Writer();
    #ifndef _USE_PARALLEL_NETCDF
        if (!mAssemble) {
            std::vector<std::vector<std::string>> allRecKeys;
            XMPI::gather(mVarNames, allRecKeys, false);
            if (XMPI::root()) {
                std::fstream fout(Parameters::sOutputDirectory + "/stations/rank_edge.txt", std::fstream::out);
                fout << "# MPI_RANK EDGE_TAG\n";
                for (int rank = 0; rank < XMPI::nproc(); rank++) {
                    for (int irec = 0; irec < allRecKeys[rank].size(); irec++) {
                        fout << rank << " " << allRecKeys[rank][irec] << "\n";
                    }
                }
            }
        }
    
        // open file on all ranks with elements
        if (numEle == 0) {
            return;
        }
        std::stringstream fname;
        fname << Parameters::sOutputDirectory + "/stations/axisem3d_surface.nc.rank" << XMPI::rank();
        mNetCDF->open(fname.str(), true);
        mNetCDF->defModeOn();
        // define time
        mNetCDF->defineVariable<double>("time_points", dimsTime);
        // define seismograms
        for (int iele = 0; iele < numEle; iele++) {
            dimsSeis[1] = nPntEdge * 3 * (mNu[iele] + 1);
            mNetCDF->defineVariable<Real>(mVarNames[iele] + "r", dimsSeis);
            mNetCDF->defineVariable<Real>(mVarNames[iele] + "i", dimsSeis);
        }
        // define theta
        mNetCDF->defineVariable<double>("theta", dimsTheta);
        // define GLL and GLJ
        mNetCDF->defineVariable<double>("GLL", dimsGLL);
        mNetCDF->defineVariable<double>("GLJ", dimsGLL);
        mNetCDF->defModeOff();
        // fill time with err values
        mNetCDF->fillConstant("time_points", dimsTime, NC_ERR_VALUE);
        // fill seismograms with err values
        for (int iele = 0; iele < numEle; iele++) {
            dimsSeis[1] = nPntEdge * 3 * (mNu[iele] + 1);
            mNetCDF->fillConstant(mVarNames[iele] + "r", dimsSeis, (Real)NC_ERR_VALUE);
            mNetCDF->fillConstant(mVarNames[iele] + "i", dimsSeis, (Real)NC_ERR_VALUE);
        }
        // fill theta
        mNetCDF->writeVariableWhole("theta", theta);
        // fill GLL and GLJ
        mNetCDF->writeVariableWhole("GLL", SpectralConstants::getP_GLL());
        mNetCDF->writeVariableWhole("GLJ", SpectralConstants::getP_GLJ());
        // source location
        mNetCDF->addAttribute("", "source_latitude", mSrcLat);
        mNetCDF->addAttribute("", "source_longitude", mSrcLon);
        mNetCDF->addAttribute("", "source_depth", mSrcDep);
        mNetCDF->addAttribute("", "source_flattening",
            Geodesy::getFlattening(Geodesy::getROuter() - mSrcDep));
        mNetCDF->addAttribute("", "surface_flattening", Geodesy::getFlattening());
        mNetCDF->addAttribute("", "radius", Geodesy::getROuter());
        mNetCDF->flush();
    #else
        // gather all variable names 
        std::vector<std::vector<std::string>> allNames;
        std::vector<std::vector<int>> allNus; 
        XMPI::gather(mVarNames, allNames, true);
        XMPI::gather(mNu, allNus, MPI_INT, true);
        
        // open file on min rank and define all variables
        std::string fname = Parameters::sOutputDirectory + "/stations/axisem3d_surface.nc";
        if (XMPI::rank() == mMinRankWithEle) {
            mNetCDF->open(fname, true);
            mNetCDF->defModeOn();
            // define time
            mNetCDF->defineVariable<double>("time_points", dimsTime);
            // define seismograms
            for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
                for (int iele = 0; iele < allNames[iproc].size(); iele++) {
                    dimsSeis[1] = nPntEdge * 3 * (allNus[iproc][iele] + 1);
                    mNetCDF->defineVariable<Real>(allNames[iproc][iele] + "r", dimsSeis);
                    mNetCDF->defineVariable<Real>(allNames[iproc][iele] + "i", dimsSeis);
                }
            }
            // define theta
            mNetCDF->defineVariable<double>("theta", dimsTheta);
            // define GLL and GLJ
            mNetCDF->defineVariable<double>("GLL", dimsGLL);
            mNetCDF->defineVariable<double>("GLJ", dimsGLL);
            mNetCDF->defModeOff();
            // fill time with err values
            mNetCDF->fillConstant("time_points", dimsTime, NC_ERR_VALUE);
            // fill seismograms with err values
            for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
                for (int iele = 0; iele < allNames[iproc].size(); iele++) {
                    dimsSeis[1] = nPntEdge * 3 * (allNus[iproc][iele] + 1);
                    mNetCDF->fillConstant<Real>(allNames[iproc][iele] + "r", dimsSeis, (Real)NC_ERR_VALUE);
                    mNetCDF->fillConstant<Real>(allNames[iproc][iele] + "i", dimsSeis, (Real)NC_ERR_VALUE);
                }
            }
            // fill theta
            mNetCDF->writeVariableWhole("theta", theta);
            // fill GLL and GLJ
            mNetCDF->writeVariableWhole("GLL", SpectralConstants::getP_GLL());
            mNetCDF->writeVariableWhole("GLJ", SpectralConstants::getP_GLJ());
            // source location
            mNetCDF->addAttribute("", "source_latitude", mSrcLat);
            mNetCDF->addAttribute("", "source_longitude", mSrcLon);
            mNetCDF->addAttribute("", "source_depth", mSrcDep);
            mNetCDF->addAttribute("", "source_flattening",
                Geodesy::getFlattening(Geodesy::getROuter() - mSrcDep));
            mNetCDF->addAttribute("", "surface_flattening", Geodesy::getFlattening());
            mNetCDF->addAttribute("", "radius", Geodesy::getROuter());
            mNetCDF->close();
        }
        XMPI::barrier();
        mNetCDF->openParallel(fname);
    #endif
    
    // record postion in nc file
    mCurrentRow = 0;
}

void SurfaceIO::finalize() {
    if (mMinRankWithEle == -1) {
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
        return;
    }
    
    // file name
    int numEle = mVarNames.size();
    std::string oneFile = Parameters::sOutputDirectory + "/stations/axisem3d_surface.nc";
    std::stringstream fname;
    fname << Parameters::sOutputDirectory + "/stations/axisem3d_surface.nc.rank" << XMPI::rank();
    std::string locFile = fname.str();
    
    // merge
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(true);
    #endif
    
    // dims
    std::vector<size_t> dimsTime;
    std::vector<size_t> dimsSeis;
    std::vector<size_t> dimsGLL;
    dimsTime.push_back(mCurrentRow);
    dimsSeis.push_back(mCurrentRow);
    dimsSeis.push_back(0);
    dimsGLL.push_back(nPntEdge);
    
    // create file 
    if (XMPI::rank() == mMinRankWithEle) {
        // read time
        NetCDF_Reader nr;
        nr.open(locFile);
        RDColX times;
        nr.read1D("time_points", times);
        
        // read theta
        RDMatXX_RM theta;
        nr.read2D("theta", theta);
        
        std::vector<size_t> dimsTheta;
        dimsTheta.push_back(theta.rows());
        dimsTheta.push_back(2);
        nr.close();
    
        // create file and write time
        NetCDF_Writer nw;
        nw.open(oneFile, true);
        nw.defModeOn();
        nw.defineVariable<double>("time_points", dimsTime);
        nw.defineVariable<double>("theta", dimsTheta);
        nw.defineVariable<double>("GLL", dimsGLL);
        nw.defineVariable<double>("GLJ", dimsGLL);
        nw.defModeOff();
        nw.writeVariableWhole("time_points", times);
        nw.writeVariableWhole("theta", theta);
        // GLL and GLJ
        nw.writeVariableWhole("GLL", SpectralConstants::getP_GLL());
        nw.writeVariableWhole("GLJ", SpectralConstants::getP_GLJ());
        // source location
        nw.addAttribute("", "source_latitude", mSrcLat);
        nw.addAttribute("", "source_longitude", mSrcLon);
        nw.addAttribute("", "source_depth", mSrcDep);
        nw.addAttribute("", "source_flattening",
            Geodesy::getFlattening(Geodesy::getROuter() - mSrcDep));
        nw.addAttribute("", "surface_flattening", Geodesy::getFlattening());
        nw.addAttribute("", "radius", Geodesy::getROuter());
        nw.close();
    }
    XMPI::barrier();
    
    // write seismograms
    for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
        if (iproc == XMPI::rank() && numEle > 0) {
            // open
            NetCDF_Reader nr;
            nr.open(locFile);
            NetCDF_Writer nw;
            nw.open(oneFile, false);
    
            // create variable
            nw.defModeOn();
            for (int iele = 0; iele < numEle; iele++) {
                dimsSeis[1] = nPntEdge * 3 * (mNu[iele] + 1);
                nw.defineVariable<Real>(mVarNames[iele] + "r", dimsSeis);
                nw.defineVariable<Real>(mVarNames[iele] + "i", dimsSeis);
            }
            nw.defModeOff();
    
            // read and write seismograms
            for (int iele = 0; iele < numEle; iele++) {
                RMatXX_RM seis_r;
                RMatXX_RM seis_i;
                nr.read2D(mVarNames[iele] + "r", seis_r);
                nr.read2D(mVarNames[iele] + "i", seis_i);
                nw.writeVariableWhole(mVarNames[iele] + "r", seis_r);
                nw.writeVariableWhole(mVarNames[iele] + "i", seis_i);
            }
    
            // close
            nw.close();
            nr.close();
        }
        XMPI::barrier();
    } 
    
    // delete local files
    if (numEle > 0) {
        std::remove(locFile.c_str());
    }
    
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(false);
    #endif
}

void SurfaceIO::dumpToFile(const std::vector<CMatXX_RM> &bufferDisp, 
    const RDColX &bufferTime, int bufferLine) {
    int numEle = mVarNames.size();
    if (bufferLine == 0) {
        return;
    }
    
    // write time
    std::vector<size_t> start;
    std::vector<size_t> count;
    start.push_back(mCurrentRow);
    count.push_back(bufferLine);
    
    // record postion in nc file
    mCurrentRow += bufferLine;
    
    #ifndef _USE_PARALLEL_NETCDF
        if (numEle == 0) {
            return;
        }  
        mNetCDF->writeVariableChunk("time_points", 
            bufferTime.topRows(bufferLine), start, count);
    #else
        if (XMPI::rank() == mMinRankWithEle) {
            mNetCDF->writeVariableChunk("time_points", 
                bufferTime.topRows(bufferLine), start, count);
        }
    #endif
    
    // write seismograms
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(true);
    #endif
    start.push_back(0);
    count.push_back(0);
    for (int iele = 0; iele < numEle; iele++) {
        count[1] = nPntEdge * 3 * (mNu[iele] + 1); 
        RMatXX_RM seis_r = bufferDisp[iele].block(0, 0, bufferLine, count[1]).real();
        RMatXX_RM seis_i = bufferDisp[iele].block(0, 0, bufferLine, count[1]).imag();
        mNetCDF->writeVariableChunk(mVarNames[iele] + "r", seis_r, start, count);
        mNetCDF->writeVariableChunk(mVarNames[iele] + "i", seis_i, start, count);
    }
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(false);
    #endif
    
    // flush to disk
    mNetCDF->flush();
}

