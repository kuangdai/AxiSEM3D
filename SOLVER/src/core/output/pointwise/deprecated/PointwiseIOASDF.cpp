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

#include <fstream>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

void PointwiseIOASDF::initialize(int totalRecordSteps, int bufferSize, 
    const std::string &components, const std::vector<PointwiseInfo> &receivers) {
    // number
    int numRec = receivers.size();
    if (numRec == 0) {
        return;
    }    
        
    // file
    mNetCDF = new NetCDF_Writer();
    
    // open file on all ranks
    std::stringstream fname;
    fname << Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.asdf.h5.rank" << XMPI::rank();
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
        mVarNames[irec] += "." + components;
        mNetCDF->defineVariable(mVarNames[irec], dims, (Real)-1.2345);
        // get info for station XML
        mNetworks.push_back(receivers[irec].mNetwork);
        mNames.push_back(receivers[irec].mName);
        mLats.push_back(receivers[irec].mLat);
        mLons.push_back(receivers[irec].mLon);
        mDeps.push_back(receivers[irec].mDep);
    }
    mComponents = components;
    
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
    std::string oneFile = Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.asdf.h5";
    std::stringstream fname;
    fname << Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.asdf.h5.rank" << XMPI::rank();
    std::string locFile = fname.str();
    
    // merge
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(true);
    #endif
    
    // find the maximum rank with receivers
    int lrank = -1;
    if (numRec > 0) {
        lrank = XMPI::rank();
    }
    int maxRank = XMPI::max(lrank);
    
    // define file
    std::string sourceID;
    std::string sourceT0_UTC;
    if (XMPI::rank() == maxRank) {
        // create file
        NetCDF_Writer nw;
        nw.open(oneFile, true);
        
        // create groups
        nw.createGroup("AuxiliaryData");
        nw.createGroup("Provenance");
        nw.createGroup("Waveforms");
        nw.addAttributeString("", "file_format", "ASDF");
        nw.addAttributeString("", "file_format_version", "1.0.0");
        
        // create QuakeML
        createQuakeML(nw, sourceID, sourceT0_UTC);
        
        // close
        nw.close();
    }
    XMPI::bcast(sourceID, maxRank);
    XMPI::bcast(sourceT0_UTC, maxRank);
    
    for (int iproc = 0; iproc < XMPI::nproc(); iproc++) {
        if (iproc == XMPI::rank() && numRec > 0) {
            // reader
            NetCDF_Reader nr;
            nr.open(locFile);
            
            // time
            double tfirst = -1.2345;
            double tlastt = -1.2345;
            double sampling_rate = -1.2345;
            RColX times;
            nr.read1D("time_points", times);
            if (times.size() > 0) {
                tfirst = times(0);
                tlastt = times(times.size() - 1);
            } 
            if (times.size() > 1) {
                sampling_rate = 1. / (times(1) - times(0));
            }
            
            // UTC
            boost::posix_time::ptime utcSource = UTCfromString(sourceT0_UTC);
            long secFirst = (long)tfirst;
            long secLastt = (long)tlastt;
            long fracFirst = (long)((tfirst - secFirst * 1.) * 1e3);
            long fracLastt = (long)((tlastt - secLastt * 1.) * 1e3);
            boost::posix_time::ptime utcFirst(utcSource + boost::posix_time::time_duration(0, 0, secFirst, fracFirst));
            boost::posix_time::ptime utcLastt(utcSource + boost::posix_time::time_duration(0, 0, secLastt, fracLastt));
            std::string tFirstStr = UTCToString(utcFirst, false);
            std::string tLasttStr = UTCToString(utcLastt, false);
            boost::posix_time::ptime utcZero = UTCfromString("1970-01-01T00:00:00");
            long long tFirstLong = (utcFirst - utcZero).total_nanoseconds();
            
            // read and write seismograms
            NetCDF_Writer nw;
            nw.open(oneFile, false);
            nw.goToGroup("Waveforms");
            for (int irec = 0; irec < numRec; irec++) {
                std::string key = mNetworks[irec] + "." + mNames[irec];
                
                // create group
                nw.createGroup(key);
                nw.goToGroup(key);
                
                // create station XML
                createStationML(nw, irec, sourceID);
                
                // read seis
                RMatXX_RM seis;
                nr.read2D(mVarNames[irec], seis);
                // write seis
                std::vector<size_t> dims;
                dims.push_back(seis.rows());
                for (int i = 0; i < 3; i++) {
                    std::string varName = key + ".." + mComponents.c_str()[i];
                    varName += "__" + tFirstStr + "__" + tLasttStr + "__synthetic";
                    nw.defineVariable(varName, dims, (Real)-1.2345);
                    nw.writeVariableWhole(varName, seis.col(i).eval());
                    // add attributes
                    nw.addAttributeString(varName, "event_id", sourceID);
                    nw.addAttribute(varName, "sampling_rate", sampling_rate);
                    nw.addAttribute(varName, "starttime", tFirstLong);
                }
                
                // back to root
                nw.goToFileRoot();
                nw.goToGroup("Waveforms");
            }
            // close writer
            nw.close();
            
            // close reader
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

void PointwiseIOASDF::createQuakeML(NetCDF_Writer &nw, std::string &sourceID, std::string &sourceT0_UTC) {
    std::string quakeStr;
    std::fstream fs(Parameters::sInputDirectory + "/ASDF/quake.xml");
    if (fs) {
        // user provided
        quakeStr = std::string((std::istreambuf_iterator<char>(fs)), 
            std::istreambuf_iterator<char>());
        
        // try to find event publicID
        std::string eventStr = "<event publicID=\"";
        size_t pos = quakeStr.find(eventStr);
        if (pos == std::string::npos) {
            throw std::runtime_error("PointwiseIOASDF::createQuakeML || "
                "Unable to find event publicID in user-provided quake.xml.");
        }
        pos += eventStr.length();
        size_t eventStart = pos;
        
        std::string eventStrEnd = "\">";
        pos = quakeStr.find(eventStrEnd, pos);
        if (pos == std::string::npos) {
            throw std::runtime_error("PointwiseIOASDF::createQuakeML || "
                "Unable to find event publicID in user-provided quake.xml.");
        }
        sourceID = quakeStr.substr(eventStart, pos - eventStart);
        pos += eventStrEnd.length();
        
        // try to find event time
        std::string originStr = "<origin publicID=\"";
        pos = quakeStr.find(originStr, pos);
        if (pos == std::string::npos) {
            throw std::runtime_error("PointwiseIOASDF::createQuakeML || "
                "Unable to find origin publicID in user-provided quake.xml.");
        }
        pos += originStr.length();
        
        std::string timeStr = "<time>";
        pos = quakeStr.find(timeStr, pos);
        if (pos == std::string::npos) {
            throw std::runtime_error("PointwiseIOASDF::createQuakeML || "
                "Unable to find \"<time>\" under \"<origion>\" in user-provided quake.xml.");
        }
        pos += timeStr.length();
        
        std::string valueStr = "<value>";
        pos = quakeStr.find(valueStr, pos);
        if (pos == std::string::npos) {
            throw std::runtime_error("PointwiseIOASDF::createQuakeML || "
                "Unable to find \"<value>\" under \"<time>\" under \"<origion>\" in user-provided quake.xml.");
        }
        pos += valueStr.length();
        size_t timeStart = pos;
        
        std::string timeStrEnd = "</value>";
        pos = quakeStr.find(timeStrEnd, pos);
        if (pos == std::string::npos) {
            throw std::runtime_error("PointwiseIOASDF::createQuakeML || "
                "Unable to find \"</value>\" under \"<time>\" under \"<origion>\" in user-provided quake.xml.");
        }
        sourceT0_UTC = quakeStr.substr(timeStart, pos - timeStart);
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
        boost::replace_first(quakeStr, "@PATH_TO_SOURCE_FILE@", mSourceFile);
        sourceID = mSourceFile;
        sourceT0_UTC = "1970-01-01T00:00:00";
    }
    fs.close();
    nw.writeStringInByte("QuakeML", quakeStr);
}

void PointwiseIOASDF::createStationML(NetCDF_Writer &nw, int irec, const std::string &sourceID) {
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
        boost::replace_first(stationStr, "@SOURCE@", boost::lexical_cast<std::string>(sourceID));
    }
    fs.close();
    nw.writeStringInByte("StationXML", stationStr);
}

boost::posix_time::ptime PointwiseIOASDF::UTCfromString(const std::string &utcStr) {
    std::string processed(utcStr);
    boost::replace_first(processed, "Z", "");
    boost::replace_first(processed, "GMT", "");
    std::vector<std::string> date_time = Parameters::splitString(processed, "T");
    // date
    std::vector<std::string> ymd = Parameters::splitString(date_time[0], "-");
    long year = boost::lexical_cast<long>(ymd[0]);
    long month = boost::lexical_cast<long>(ymd[1]);
    long day = boost::lexical_cast<long>(ymd[2]);
    // time
    std::vector<std::string> hms = Parameters::splitString(date_time[1], ":");
    long hour = boost::lexical_cast<long>(hms[0]);
    long min = boost::lexical_cast<long>(hms[1]);
    std::vector<std::string> sec_frac = Parameters::splitString(hms[2], ".");
    long second = boost::lexical_cast<long>(sec_frac[0]);
    long frac = 0;
    if (sec_frac.size() > 1) {
        if (sec_frac[1].length() > 0) {
            frac = boost::lexical_cast<long>(sec_frac[1]);
        }
    }
    return boost::posix_time::ptime(boost::gregorian::date(year, month, day), 
        boost::posix_time::time_duration(hour, min, second, frac));
}

std::string PointwiseIOASDF::UTCToString(const boost::posix_time::ptime &utc, bool printfrac) {
    // date
    boost::gregorian::date myDate = utc.date();
    long year = myDate.year();
    long month = myDate.month();
    long day = myDate.day();
    // time
    boost::posix_time::time_duration myTime = utc.time_of_day();
    long hour = myTime.hours();
    long min = myTime.minutes();
    long second = myTime.seconds();
    long frac = myTime.fractional_seconds();
    std::stringstream ss;
    ss << std::setfill ('0') << std::setw(4) << year << "-";
    ss << std::setfill ('0') << std::setw(2) << month << "-";
    ss << std::setfill ('0') << std::setw(2) << day << "T";
    ss << std::setfill ('0') << std::setw(2) << hour << ":";
    ss << std::setfill ('0') << std::setw(2) << min << ":";
    ss << std::setfill ('0') << std::setw(2) << second;
    if (frac > 0 && printfrac) {
        ss << "." << frac;
    }
    return ss.str();
}

