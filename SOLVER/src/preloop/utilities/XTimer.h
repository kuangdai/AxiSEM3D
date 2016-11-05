// XTimer.h
// created by Kuangdai on 5-Nov-2016 
// Recursive timer

#pragma once

#include <boost/timer/timer.hpp>
#include <vector>
#include <fstream>

class XTimer {
public:
    static void initialize(std::string fileName, int nLevels);
    static void finalize();
    static void begin(std::string procName, int level);
    static void end(std::string procName, int level);
    static void enable() {mEnabled = true;};
    static void disable() {mEnabled = false;};
    
private:
    static std::vector<boost::timer::cpu_timer> mTimers;
    static std::fstream mFile; 
    static bool mEnabled;
};

