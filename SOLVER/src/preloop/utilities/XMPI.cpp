// XMPI.cpp
// created by Kuangdai on 28-Jun-2016 
// mpi interfaces

#include "XMPI.h"
#include <iomanip>
#include "Parameters.h"

XMPI::root_cout XMPI::cout;
std::string XMPI::endl = "\n";

void XMPI::initialize(int argc, char *argv[]) {
    #ifndef _SERIAL_BUILD
        MPI_Init(&argc, &argv);
    #endif
    
    // find path of executable
    std::string argv0(argv[0]);
    std::size_t found = argv0.find_last_of("/\\");
    if (found == std::string::npos) {
        argv0 = "./" + argv0;
        found = argv0.find_last_of("/\\");
    }
    std::string execDirectory = argv0.substr(0, found);
    
    // so far, this problem happens with valgrind only 
    if (execDirectory.length() == 0) {
        execDirectory = ".";
    }
    
    // input and output
    Parameters::sInputDirectory = execDirectory + "/input";
    Parameters::sOutputDirectory = execDirectory + "/output";
    if (XMPI::root()) {
        if (!dirExists(Parameters::sInputDirectory)) {
            throw std::runtime_error("XMPI::initialize || Missing input directory: ||" 
                + Parameters::sInputDirectory);
        }
        mkdir(Parameters::sOutputDirectory);
        mkdir(Parameters::sOutputDirectory + "/stations");    
        mkdir(Parameters::sOutputDirectory + "/plots");
        mkdir(Parameters::sOutputDirectory + "/develop");
    }
}

void XMPI::finalize() {
    #ifndef _SERIAL_BUILD
        MPI_Finalize();
    #endif
}

void XMPI::printException(const std::exception &e) {
    std::string head = " AXISEM3D ABORTED UPON RUNTIME EXCEPTION ";
    std::string what = e.what();
    std::vector<std::string> strs = Parameters::splitString(what, "|");
    std::string src, msg;
    int nmsg = 0;
    if (strs.size() >= 2) {
        src = "FROM: " + boost::trim_copy(strs[0]);
        nmsg = src.length();
        msg = "WHAT: ";
        for (int i = 1; i < strs.size(); i++) {
            msg += boost::trim_copy(strs[i]);
            if (i != strs.size() - 1) {
                msg += "\n      ";
            }
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
    if (msg.length() > 0) {
        XMPI::cout << msg << XMPI::endl;
    }
    XMPI::cout << std::setfill('*') << std::setw(nstar * 2 + head.length()) << "";
    XMPI::cout << XMPI::endl << XMPI::endl;
}

void XMPI::bcast(int *buffer, int size, int src) {
    #ifndef _SERIAL_BUILD
        MPI_Bcast(buffer, size, MPI_INT, src, MPI_COMM_WORLD);
    #endif
}

void XMPI::bcast(double *buffer, int size, int src) {
    #ifndef _SERIAL_BUILD
        MPI_Bcast(buffer, size, MPI_DOUBLE, src, MPI_COMM_WORLD);
    #endif
}

void XMPI::bcast(float *buffer, int size, int src) {
    #ifndef _SERIAL_BUILD
        MPI_Bcast(buffer, size, MPI_FLOAT, src, MPI_COMM_WORLD);
    #endif
}

void XMPI::bcast(std::complex<double> *buffer, int size, int src) {
    #ifndef _SERIAL_BUILD
        MPI_Bcast(buffer, size, MPI_C_DOUBLE_COMPLEX, src, MPI_COMM_WORLD);
    #endif
}

void XMPI::bcast(std::complex<float> *buffer, int size, int src) {
    #ifndef _SERIAL_BUILD
        MPI_Bcast(buffer, size, MPI_C_FLOAT_COMPLEX, src, MPI_COMM_WORLD);
    #endif
}

void XMPI::bcast(char *buffer, int size, int src) {
    #ifndef _SERIAL_BUILD
        MPI_Bcast(buffer, size, MPI_CHAR, src, MPI_COMM_WORLD);
    #endif
}

void XMPI::bcast(int &buffer, int src) {
    #ifndef _SERIAL_BUILD
        MPI_Bcast(&buffer, 1, MPI_INT, src, MPI_COMM_WORLD);
    #endif
}

void XMPI::bcast(double &buffer, int src) {
    #ifndef _SERIAL_BUILD
        MPI_Bcast(&buffer, 1, MPI_DOUBLE, src, MPI_COMM_WORLD);
    #endif
}

void XMPI::bcast(float &buffer, int src) {
    #ifndef _SERIAL_BUILD
        MPI_Bcast(&buffer, 1, MPI_FLOAT, src, MPI_COMM_WORLD);
    #endif
}

void XMPI::bcast(std::string &str, int src) {
    #ifndef _SERIAL_BUILD
        int size = 0;
        if (rank() == src) {
            size = str.size();
        }
        bcast(size, src);
        char *cstr = new char[size + 1];
        if (rank() == src) {
            const char *ostr = str.c_str();
            for (int i = 0; i < size; i++) {
                cstr[i] = ostr[i];
            }
        } 
        bcast(cstr, size, src);
        if (rank() != src) {
            cstr[size] = 0;
            str = std::string(cstr);
        }
        delete [] cstr;
    #endif
}

void XMPI::bcast(std::vector<std::string> &buffer, int src) {
    #ifndef _SERIAL_BUILD
        // length
        int num = 0;
        int maxLen = -1;
        if (rank() == src) {
            num = buffer.size();
            for (int i = 0; i < num; i++) {
                maxLen = std::max(maxLen, (int)buffer[i].size());
            }
            maxLen += 1;
        }
        bcast(num, src);
        bcast(maxLen, src);
        
        // bcast data
        char *all_str = new char[num * maxLen];
        for (int i = 0; i < num * maxLen; i++) {
            all_str[i] = 0;
        }
        if (rank() == src) {
            for (int i = 0; i < num; i++) {
                for (int j = 0; j < buffer[i].size(); j++) {
                    all_str[i * maxLen + j] = buffer[i][j];
                }
            }
        }
        bcast(all_str, num * maxLen, src);
        
        // cast back to array
        if (rank() != src) {
            buffer.clear();
            for (int i = 0; i < num; i++) {
                std::string var(&all_str[i * maxLen]);
                buffer.push_back(var);
            }
        }
        
        delete [] all_str;
    #endif
}

void XMPI::min(const std::vector<int> &value, std::vector<int> &minimum) {
    #ifndef _SERIAL_BUILD
        MPI_Allreduce(value.data(), minimum.data(), value.size(), MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    #endif
}

int XMPI::min(const int &value) {
    #ifndef _SERIAL_BUILD
        int minimum;
        MPI_Allreduce(&value, &minimum, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        return minimum;
    #else
        return value;
    #endif
}

double XMPI::min(const double &value) {
    #ifndef _SERIAL_BUILD
        double minimum;
        MPI_Allreduce(&value, &minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        return minimum;
    #else
        return value;
    #endif
}

int XMPI::max(const int &value) {
    #ifndef _SERIAL_BUILD
        int minimum;
        MPI_Allreduce(&value, &minimum, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        return minimum;
    #else
        return value;
    #endif
}

double XMPI::max(const double &value) {
    #ifndef _SERIAL_BUILD
        double minimum;
        MPI_Allreduce(&value, &minimum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return minimum;
    #else
        return value;
    #endif
}

int XMPI::sum(const int &value) {
    #ifndef _SERIAL_BUILD
        int minimum;
        MPI_Allreduce(&value, &minimum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        return minimum;
    #else
        return value;
    #endif
}

double XMPI::sum(const double &value) {
    #ifndef _SERIAL_BUILD
        double minimum;
        MPI_Allreduce(&value, &minimum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return minimum;
    #else
        return value;
    #endif
}

void XMPI::sumVector(std::vector<double> &value) {
    #ifndef _SERIAL_BUILD
        std::vector<double> total(value.size());
        MPI_Allreduce(value.data(), total.data(), value.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        value = total;
    #endif
}

void XMPI::gather(int buf, std::vector<int> &all_buf, bool all) {
    #ifndef _SERIAL_BUILD
        if (all) {
            all_buf.resize(XMPI::nproc());
            MPI_Allgather(&buf, 1, MPI_INT, all_buf.data(), 1, MPI_INT, MPI_COMM_WORLD);
        } else {
            if (root()) {
                all_buf.resize(XMPI::nproc());
            }
            MPI_Gather(&buf, 1, MPI_INT, all_buf.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
    #else
        all_buf.clear();
        all_buf.push_back(buf);
    #endif
}

void XMPI::gather(const std::string &buf, std::vector<std::string> &all_buf, bool all) {
    #ifndef _SERIAL_BUILD
        // size
        int size = buf.size() + 1;
        std::vector<int> all_size;
        gather(size, all_size, all);
        int total_size = 0;
        for (auto &n : all_size) {
            total_size += n;
        }
        
        // displacement: where to insert to the global array
        int nproc = XMPI::nproc();
        std::vector<int> disp(nproc, 0);
        if (all || root()) {
            for (int i = 1; i < nproc; i++) {
                for (int j = 0; j < i; j++) {
                    disp[i] += all_size[j];
                }
            }
        }

        char *all_cstr;
        if (all || root()) {
            all_cstr = new char[total_size];
        }
        if (all) {
            MPI_Allgatherv(buf.c_str(), size, MPI_CHAR, all_cstr, all_size.data(), disp.data(), MPI_CHAR, MPI_COMM_WORLD);
        } else {
            MPI_Gatherv(buf.c_str(), size, MPI_CHAR, all_cstr, all_size.data(), disp.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
        }
        if (all || root()) {
            all_buf.clear();
            int pos = 0;
            for (int i = 0; i < nproc; i++) {
                all_buf.push_back(std::string(&all_cstr[pos]));
                pos += all_size[i];
            }
            delete [] all_cstr;
        }
        
    #else
        all_buf.clear();
        all_buf.push_back(buf);
    #endif
}

void XMPI::gather(const std::vector<std::string> &buf, 
    std::vector<std::vector<std::string>> &all_buf, 
    bool all) {
    
    std::string loc_buf_all = "";
    std::vector<int> loc_buf_size(buf.size(), 0);
    for (int i = 0; i < buf.size(); i++) {
        loc_buf_all += buf[i];
        loc_buf_size[i] = buf[i].size();
    }
    
    std::vector<std::string> all_buf_flat;
    std::vector<std::vector<int>> all_size;
    gather(loc_buf_all, all_buf_flat, all);
    gather(loc_buf_size, all_size, MPI_INT, all);
    
    if (all || root()) {
        all_buf.clear();
        for (int i = 0; i < all_size.size(); i++) {
            std::vector<std::string> str;
            int pos = 0;
            for (int j = 0; j < all_size[i].size(); j++) {
                int length = all_size[i][j];
                str.push_back(all_buf_flat[i].substr(pos, length));
                pos += length;
            }
            all_buf.push_back(str);
        }
    }
}


/////////////////
// file system 
/////////////////

extern "C" {
    #include <sys/types.h>
    #include <sys/stat.h>
};

bool XMPI::dirExists(const std::string &path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        return false;
    } else if (info.st_mode & S_IFDIR) {
        return true;
    } else {
        return false;
    }
}

void XMPI::mkdir(const std::string &path) {
    if (!dirExists(path)) {
        ::mkdir(path.c_str(), ACCESSPERMS);
    }
}
