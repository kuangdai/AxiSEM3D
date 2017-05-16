// MultilevelTimer.h
// created by Kuangdai on 5-Nov-2016 
// multi-level timer

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


class MultilevelTimer {
public:
    static void initialize(const std::string &fileName, int nLevels);
    static void finalize();
    static void begin(const std::string &procName, int level, bool barrier = false);
    static void end(const std::string &procName, int level, bool barrier = true);
    static void pause(int level);
    static void resume(int level);
    
    static void enable() {openFile(); mEnabled = true;};
    static void disable() {mEnabled = false;};
    
private:
    static void openFile();
    static std::string mFileName;
    static std::fstream mFile; 
    static std::vector<MyBoostTimer> mTimers;
    static bool mEnabled;
};

