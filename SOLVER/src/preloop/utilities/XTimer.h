// XTimer.h
// created by Kuangdai on 5-Nov-2016 
// Recursive timer

#pragma once

#include <vector>
#include <fstream>
#include <chrono>

class MyBoostTimer {
public:
    void stop();
    void resume();
    void clear();
    double elapsed(); // in seconds
    void start() {clear(); resume();};
    
    bool running() {return mTimePoints.size() % 2 == 1;};
    static double getClockResolution();
    
private:
    std::vector<std::chrono::high_resolution_clock::time_point> mTimePoints;
};


class XTimer {
public:
    static void initialize(const std::string &fileName, size_t nLevels);
    static void finalize();
    static void begin(const std::string &procName, size_t level, bool barrier = false);
    static void end(const std::string &procName, size_t level, bool barrier = true);
    static void pause(size_t level);
    static void resume(size_t level);
    
    static void enable() {openFile(); mEnabled = true;};
    static void disable() {mEnabled = false;};
    
private:
    static void openFile();
    static std::string mFileName;
    static std::fstream mFile; 
    static std::vector<MyBoostTimer> mTimers;
    static bool mEnabled;
};

