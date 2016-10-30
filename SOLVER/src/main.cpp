// main.cpp
// created by Kuangdai on 26-Mar-2016 
// main 

#include "axisem.h"

extern "C" void set_ftz();

int main(int argc, char *argv[]) {
    
    // denormal float handling 
    set_ftz();
    
    // axisem main
    return axisem_main(argc, argv);
    
}
