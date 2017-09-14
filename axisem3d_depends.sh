#!/bin/bash
# axisem3d_depends.sh
# created by Kuangdai on 30-Oct-2016 
# one-click installation of AxiSEM3D dependencies


############################## edit within this box ##############################
# installation path of conda
export MY_CONDA_INSTALL_DIR=$HOME/anaconda

# your bash_profile
# Mac OS X default
export MY_BASH_PROFILE_OSX=$HOME/.bash_profile
# Ubuntu default
export MY_BASH_PROFILE_UBT=$HOME/.bashrc

# What do you have already? And where are they?
# boost
export MY_BOOST_READY=false
export MY_BOOST_DIR=""
# eigen3 (make sure it is above 3.3-rc1, very important)
export MY_EIGEN3_READY=false
export MY_EIGEN3_DIR=""
# fftw (make sure it contains both double- and single-precision)
export MY_FFTW_READY=false
export MY_FFTW_DIR=""
# Metis (make sure it is built with 32-bit)
export MY_METIS_READY=false
export MY_METIS_DIR=""
# NetCDF
export MY_NETCDF_READY=false
export MY_NETCDF_DIR=""

# where you would like to download boost and eigen3
export BOOST_INSTALL_DIR=$HOME/axisem3d_depends/boost 
export EIGEN3_INSTALL_DIR=$HOME/axisem3d_depends/eigen3 
############################## edit within this box ##############################

##### Boost #####
if [ $MY_BOOST_READY == false ]; then
    mkdir -p $BOOST_INSTALL_DIR
    if [ ! -f $BOOST_INSTALL_DIR/boost_1_63_0.tar.bz2 ]; then
        wget -P $BOOST_INSTALL_DIR https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.bz2
    fi
    tar -jxf $BOOST_INSTALL_DIR/boost_1_63_0.tar.bz2 --strip 1 -C $BOOST_INSTALL_DIR
    export MY_BOOST_DIR=$BOOST_INSTALL_DIR
    export MY_BOOST_READY=true
fi

##### EIGEN3 #####
if [ $MY_EIGEN3_READY == false ]; then
    mkdir -p $EIGEN3_INSTALL_DIR
    if [ ! -f $EIGEN3_INSTALL_DIR/3.3-rc1.tar.bz2 ]; then
        wget -P $EIGEN3_INSTALL_DIR http://bitbucket.org/eigen/eigen/get/3.3-rc1.tar.bz2
    fi
    tar -jxf $EIGEN3_INSTALL_DIR/3.3-rc1.tar.bz2 --strip 1 -C $EIGEN3_INSTALL_DIR
    export MY_EIGEN3_DIR=$EIGEN3_INSTALL_DIR
    export MY_EIGEN3_READY=true
fi

##### FFTW #####
if [ $MY_FFTW_READY == false ]; then
    conda install -c menpo fftw=3.3.4
    export MY_FFTW_DIR=$MY_CONDA_INSTALL_DIR
    export MY_FFTW_READY=true
fi

##### METIS #####
if [ $MY_METIS_READY == false ]; then
    conda install -c menpo metis=5.1.0
    export MY_METIS_DIR=$MY_CONDA_INSTALL_DIR
    export MY_METIS_READY=true
fi

##### NetCDF #####
if [ $MY_NETCDF_READY == false ]; then
    conda install netcdf4
    export MY_NETCDF_DIR=$MY_CONDA_INSTALL_DIR
    export MY_NETCDF_READY=true
fi

############### add ROOTs for cmake of AxiSEM3D ###############
export MY_AXISEM_ROOTS=$HOME/.axisem3d_roots
if [ -f $MY_AXISEM_ROOTS ]; then
    rm $MY_AXISEM_ROOTS
fi
echo "export BOOST_ROOT=$MY_BOOST_DIR" >> $MY_AXISEM_ROOTS
echo "export EIGEN3_ROOT=$MY_EIGEN3_DIR" >> $MY_AXISEM_ROOTS
echo "export FFTW_ROOT=$MY_FFTW_DIR" >> $MY_AXISEM_ROOTS
echo "export METIS_ROOT=$MY_METIS_DIR" >> $MY_AXISEM_ROOTS
echo "export NETCDF_ROOT=$MY_NETCDF_DIR" >> $MY_AXISEM_ROOTS
. $MY_AXISEM_ROOTS

# add to bash profile
if [ -f $MY_BASH_PROFILE_OSX ]; then
    if ! grep -Fxq ". $MY_AXISEM_ROOTS" $MY_BASH_PROFILE_OSX 
    then
        echo ". $MY_AXISEM_ROOTS" >> $MY_BASH_PROFILE_OSX
    fi
fi

if [ -f $MY_BASH_PROFILE_UBT ]; then
    if ! grep -Fxq ". $MY_AXISEM_ROOTS" $MY_BASH_PROFILE_UBT 
    then
        echo ". $MY_AXISEM_ROOTS" >> $MY_BASH_PROFILE_UBT
    fi
fi

