// XMPI.cpp
// created by Kuangdai on 28-Jun-2016 
// mpi interfaces

#include "XMPI.h"
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include "Parameters.h"

#ifndef _SERIAL_BUILD
    boost::mpi::environment *XMPI::sEnv = 0;
    boost::mpi::communicator *XMPI::sWorld = 0;
#endif

XMPI::root_cout XMPI::cout;
std::string XMPI::endl = "\n";

void XMPI::initialize(int argc, char *argv[]) {
    #ifndef _SERIAL_BUILD
        XMPI::sEnv = new boost::mpi::environment(argc, argv);
        XMPI::sWorld = new boost::mpi::communicator();
    #endif
    std::string argv0(argv[0]);
    std::string execDirectory = argv0.substr(0, argv0.length() - 9);
    if (boost::algorithm::ends_with(execDirectory, "/.")) 
        execDirectory = execDirectory.substr(0, execDirectory.length() - 2);
    // when launched by some debuggers such as valgrind, boost cannot find exe directory
    if (execDirectory.length() == 0) execDirectory = ".";
    Parameters::sInputDirectory = execDirectory + "/input";
    Parameters::sOutputDirectory = execDirectory + "/output";
    if (!boost::filesystem::exists(Parameters::sInputDirectory)) 
        throw std::runtime_error("XMPI::initialize || Missing input directory: ||" + Parameters::sInputDirectory);
    if (!boost::filesystem::exists(Parameters::sOutputDirectory)) 
        boost::filesystem::create_directory(Parameters::sOutputDirectory);
}

void XMPI::finalize() {
    #ifndef _SERIAL_BUILD
        delete sWorld;
        delete sEnv;
    #endif
}

void XMPI::printException(const std::exception &e) {
    std::string head = " AXISEM3D ABORTED UPON RUNTIME EXCEPTION ";
    std::string what = e.what();
    std::vector<std::string> strs;
    boost::trim_if(what, boost::is_any_of("\t "));
    boost::split(strs, what, boost::is_any_of("|"), boost::token_compress_on);
    std::string src, msg;
    int nmsg = 0;
    if (strs.size() >= 2) {
        src = "FROM: " + boost::trim_copy(strs[0]);
        nmsg = src.length();
        msg = "WHAT: ";
        for (int i = 1; i < strs.size(); i++) {
            msg += boost::trim_copy(strs[i]);
            if (i != strs.size() - 1) msg += "\n      ";
            nmsg = std::max(nmsg, (int)(boost::trim_copy(strs[i]).length() + 6));
        }
    } else {
        src = what;
        msg = "";
        nmsg = src.length();
    }
    nmsg += 2;
    int nstar = std::max(5, (int)(nmsg - head.length()) / 2);
    XMPI::cout << XMPI::endl << std::setfill('*') << std::setw(nstar) << "";
    XMPI::cout << head << std::setfill('*') << std::setw(nstar) << "" << XMPI::endl;
    XMPI::cout << src << XMPI::endl;
    if (msg.length() > 0) XMPI::cout << msg << XMPI::endl;
    XMPI::cout << std::setfill('*') << std::setw(nstar * 2 + head.length()) << "";
    XMPI::cout << XMPI::endl << XMPI::endl;
}

