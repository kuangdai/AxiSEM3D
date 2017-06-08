// Parameters.h
// created by Kuangdai on 28-Jun-2016 
// simulation parameters

#pragma once

#include <vector>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <memory>
#include <typeinfo>

class Parameters {
public:
    
    // register, read, broadcast
    void initReadBcast();
    
    // parse a line from input file
    void parseLine(const std::string &line);
    
    // cast and get value 
    template<typename parType>
    parType getValue(const std::string &key, int index = 0) const {
        try {
            // special care of bool
            std::string val = mKeyValues.at(key).at(index);
            if (typeid(parType) == typeid(bool)) {
                boost::to_upper<std::string>(val);
                if (val == "TRUE" || val == "YES" || val == "ON") val = "1";
                if (val == "FALSE" || val == "NO" || val == "OFF") val = "0";
            }
            return boost::lexical_cast<parType>(val);
        } catch(std::exception) {
            throw std::runtime_error("Parameters::getValue || "
                "Invalid parameter, keyword = " + key + ".");
        }
    };
    
    // get size of values
    int getSize(const std::string &key) const {
        try {
            return mKeyValues.at(key).size();
        } catch(std::exception) {
            throw std::runtime_error("Parameters::getSize || "
                "Invalid parameter, keyword = " + key + ".");
        }
    };
    
    // verbose
    std::string verbose() const;
    
    // build from input parameters
    static void buildInparam(Parameters *&par, int &verbose);
        
    // string cast tools 
    template<typename parType>
    static void castValue(parType &result, const std::string &val_in, const std::string &source) {
        try {
            // special care of bool
            std::string val = val_in;
            if (typeid(parType) == typeid(bool)) {
                boost::to_upper<std::string>(val);
                if (val == "TRUE" || val == "YES" || val == "ON") val = "1";
                if (val == "FALSE" || val == "NO" || val == "OFF") val = "0";
            }
            result = boost::lexical_cast<parType>(val);
        } catch(std::exception) {
            throw std::runtime_error("Parameters::castValue || "
                "Invalid argument encountered in " + source + ", arg = " + val_in + ".");
        }
    };
    
    static std::vector<std::string> splitString(const std::string &in, const std::string &sep);
    
    // input and output dir 
    static std::string sInputDirectory;
    static std::string sOutputDirectory;
    
private:    
    void registerAll();
    void registerPar(const std::string &key);
    void readParFile(const std::string &fname);
    
    // map of parameters 
    std::map<std::string, std::vector<std::string>> mKeyValues; 
};

