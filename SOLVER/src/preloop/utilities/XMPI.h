// XMPI.h
// created by Kuangdai on 28-Jun-2016 
// mpi interfaces

#pragma once 
#include <iostream>
#include "eigenc.h"

#ifndef _SERIAL_BUILD
    #include <boost/mpi.hpp>
    #include <boost/serialization/string.hpp>
    #include <boost/serialization/array.hpp>
    #include <boost/serialization/vector.hpp>
    #include <boost/serialization/map.hpp>
    #include <boost/serialization/complex.hpp>
#endif

class XMPI {
public:
    // initialize and finalize
    static void initialize(int argc, char *argv[]);
    static void finalize();
    
    // print exception
    static void printException(const std::exception &e);
    
    // properties
    static int nproc() {
        #ifndef _SERIAL_BUILD
            return sWorld->size();
        #else
            return 1;
        #endif
    };
    
    static int rank() {
        #ifndef _SERIAL_BUILD
            return sWorld->rank();
        #else
            return 0;
        #endif
    };
    
    static bool root() {return rank() == 0;};
    
    // barrier
    static void barrier() {
        #ifndef _SERIAL_BUILD
            sWorld->barrier();
        #endif
    };
    
    // abort
    static void abort(int err = 0) {
        #ifndef _SERIAL_BUILD
            sEnv->abort(err);
        #else
            exit(0);
        #endif
    };
    
    ////////////////////////////// broadcast //////////////////////////////
    // array
    template<typename Type>
    static void bcast(Type *buffer, int size) {
        #ifndef _SERIAL_BUILD
            boost::mpi::broadcast(*sWorld, buffer, size, 0);
        #endif
    };
    
    // single
    template<typename Type>
    static void bcast(Type &buffer) {
        #ifndef _SERIAL_BUILD
            boost::mpi::broadcast(*sWorld, buffer, 0);
        #endif
    };
    
    // single from iproc
    template<typename Type>
    static void bcastFromProc(Type &buffer, int iproc) {
        #ifndef _SERIAL_BUILD
            boost::mpi::broadcast(*sWorld, buffer, iproc);
        #endif
    };
    
    // Eigen::Matrix
    template<typename Type>
    static void bcastEigen(Type &buffer) {
        #ifndef _SERIAL_BUILD
            int row, col;
            if (root()) {row = buffer.rows(); col = buffer.cols();}
            bcast(row);
            bcast(col);
            if (!root()) buffer = Type::Zero(row, col);
            bcast(buffer.data(), buffer.size());
        #endif
    };
    
    ////////////////////////////// isend/irecv ////////////////////////////// 
    // request typdef
    #ifndef _SERIAL_BUILD
        typedef boost::mpi::request Request; 
    #else 
        typedef int Request;
    #endif
    
    // isend, only for Eigen::Matrix
    template<typename EigenMat>
    static Request isend(int dest, const EigenMat &buffer) {
        #ifndef _SERIAL_BUILD
            return sWorld->isend(dest, dest, buffer.data(), buffer.size());
        #else
            return 0;
        #endif
    };
    
    // irecv, only for Eigen::Matrix
    template<typename EigenMat>
    static Request irecv(int source, EigenMat &buffer) {
        #ifndef _SERIAL_BUILD
            return sWorld->irecv(source, rank(), buffer.data(), buffer.size());
        #else
            return 0;
        #endif
    };
    
    // wait_all
    template<typename It>
    static void wait_all(It first, It last) {
        #ifndef _SERIAL_BUILD
            boost::mpi::wait_all(first, last);
        #endif
    };
    
    ////////////////////////////// reduce //////////////////////////////
    // minimum
    template<typename Type>
    static Type min(const Type &value) {
        #ifndef _SERIAL_BUILD
            Type minimum;
            boost::mpi::all_reduce(*sWorld, value, minimum, boost::mpi::minimum<Type>());
            return minimum;
        #else
            return value;
        #endif
    };
    
    // sum
    template<typename Type>
    static Type sum(const Type &value) {
        #ifndef _SERIAL_BUILD
            Type total;
            boost::mpi::all_reduce(*sWorld, value, total, std::plus<Type>());
            return total;
        #else 
            return value;
        #endif
    };
    
    // sum Eigen::Matrix
    template<typename Type>
    static void sumEigen(Type &value) {
        #ifndef _SERIAL_BUILD
            Type total(value);
            boost::mpi::all_reduce(*sWorld, value.data(), value.size(), total.data(), 
                std::plus<typename Type::Scalar>());
            value = total;
        #endif
    };
    
    ////////////////////////////// gather //////////////////////////////
    template<typename Type>
    static std::vector<Type> all_gather(const Type &value) {
        #ifndef _SERIAL_BUILD
            std::vector<Type> total;
            boost::mpi::all_gather(*sWorld, value, total);
            return total;
        #else 
            return std::vector<Type>(value, 1);
        #endif
    };
    
    template<typename Type>
    static std::vector<Type> gather(const Type &value) {
        #ifndef _SERIAL_BUILD
            std::vector<Type> total;
            boost::mpi::gather(*sWorld, value, total, 0);
            return total;
        #else 
            return std::vector<Type>(value, 1);
        #endif
    };
        
    ////////////////////////////// stream on root //////////////////////////////
    struct root_cout {
        template<typename Type>
        root_cout &operator<<(const Type &val) {
            if (rank() == mProc) std::cout << val;
            return *this;
        };
        void resetp() {mProc = 0;};
        void setp(int proc) {
            if (proc < nproc()) mProc = proc;
            else resetp();
        };
        int mProc = 0;
    };
    static root_cout cout;
    static std::string endl;
        
private:
    #ifndef _SERIAL_BUILD
        static boost::mpi::environment *sEnv;
        static boost::mpi::communicator *sWorld;
    #endif
};

// message info
struct MessagingInfo {
    // number of procs to communicate with
    int mNProcComm;
    // procs to communicate with
    std::vector<int> mIProcComm;
    // number of local points to be communicated for each proc
    std::vector<int> mNLocalPoints;
    // indecies of local points to be communicated for each proc
    std::vector<std::vector<int>> mILocalPoints;
    // mpi requests
    std::vector<XMPI::Request> mReqSend;
    std::vector<XMPI::Request> mReqRecv;
};

// message buffer for solver
struct MessagingBuffer {
    std::vector<CColX> mBufferSend;
    std::vector<CColX> mBufferRecv;
};


