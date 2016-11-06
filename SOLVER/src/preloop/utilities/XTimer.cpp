// XTimer.cpp
// created by Kuangdai on 5-Nov-2016 
// Recursive timer

#include "XTimer.h"
#include "XMPI.h"

std::vector<boost::timer::cpu_timer> XTimer::mTimers;
std::fstream XTimer::mFile;
bool XTimer::mEnabled = false;

void XTimer::initialize(std::string fileName, int nLevels) {
    if (!XMPI::root()) return;
    mFile = std::fstream(fileName, std::fstream::out);
    mTimers.clear();
    for (int i = 0; i < nLevels; i++) {
        mTimers.push_back(boost::timer::cpu_timer());
        mTimers[i].stop();
        mTimers[i].elapsed().clear();
    }
}

void XTimer::finalize() {
    if (!XMPI::root()) return;
    mFile.close();
}

void XTimer::begin(std::string name, int level, bool barrier) {
    if (!mEnabled) return;
    if (barrier) XMPI::barrier();
    if (!XMPI::root()) return;
    for (int i = 0; i < level; i++) mFile << "    ";
    mFile << name << " begins..." << std::endl;
    mTimers[level].elapsed().clear();
    mTimers[level].start();
}

void XTimer::end(std::string name, int level, bool barrier) {
    if (!mEnabled) return;
    if (barrier) XMPI::barrier();
    if (!XMPI::root()) return;
    mTimers[level].stop();
    for (int i = 0; i < level; i++) mFile << "    ";
    mFile << name << " finishes. Elapsed seconds = " << mTimers[level].elapsed().wall / 1e9 << std::endl;
}

void XTimer::pause(int level) {
    if (!XMPI::root() || !mEnabled) return;
    mTimers[level].stop();
}

void XTimer::resume(int level) {
    if (!XMPI::root() || !mEnabled) return;
    mTimers[level].resume();
}

