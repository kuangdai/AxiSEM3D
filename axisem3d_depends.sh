#!/bin/bash
# axisem3d_depends.sh
# created by Kuangdai on 30-Oct-2016 
# one-click installation of AxiSEM3D dependencies


############################## edit within this box ##############################
# installation path of conda
export MY_CONDA_INSTALL_DIR=$HOME/anaconda

# your bash_profile
# Mac OS X default
export MY_BASH_PROFILE=$HOME/.bash_profile
# Ubuntu default
# export MY_BASH_PROFILE=$HOME/.bashrc

# What do you have already? And where are they?
# boost
export MY_BOOST_READY=false
export MY_BOOST_DIR=""
# fftw (make sure it contains both double- and single-precision)
export MY_FFTW_READY=false
export MY_FFTW_DIR=""
# Metis (make sure it is built with 32-bit)
export MY_METIS_READY=false
export MY_METIS_DIR=""
# eigen3 (make sure it is above 3.3-rc1, very important)
export MY_EIGEN3_READY=false
export MY_EIGEN3_DIR=""
# HDF5
export MY_HDF5_READY=false
export MY_HDF5_DIR=""

# where you would like to install eigen3
export EIGEN3_INSTALL_DIR=$HOME/axisem3d_depends/eigen3 
############################## edit within this box ##############################

##### Boost #####
if [ $MY_BOOST_READY == false ]; then
    conda install -c anaconda boost=1.61.0
    export MY_BOOST_DIR=$MY_CONDA_INSTALL_DIR
    export MY_BOOST_READY=true
fi

##### FFTW #####
if [ $MY_FFTW_READY == false ]; then
    conda install -c salilab fftw=3.3.4
    export MY_FFTW_DIR=$MY_CONDA_INSTALL_DIR
    export MY_FFTW_READY=true
fi

##### METIS #####
if [ $MY_METIS_READY == false ]; then
    conda install -c menpo metis=5.1.0
    export MY_METIS_DIR=$MY_CONDA_INSTALL_DIR
    export MY_METIS_READY=true
fi

##### EIGEN3 #####
if [ $MY_EIGEN3_READY == false ]; then
    mkdir -p $EIGEN3_INSTALL_DIR
    if [ ! -f $EIGEN3_INSTALL_DIR/3.3-rc1.tar.bz2 ]; then
        wget -P $EIGEN3_INSTALL_DIR http://bitbucket.org/eigen/eigen/get/3.3-rc1.tar.bz2
    fi
    tar -jxvf $EIGEN3_INSTALL_DIR/3.3-rc1.tar.bz2 --strip 1 -C $EIGEN3_INSTALL_DIR
    export MY_EIGEN3_DIR=$EIGEN3_INSTALL_DIR
    export MY_EIGEN3_READY=true
fi

##### HDF5 #####
if [ $MY_HDF5_READY == false ]; then
    conda install -c conda-forge hdf5=1.8.17
    export MY_HDF5_DIR=$MY_CONDA_INSTALL_DIR
    export MY_HDF5_READY=true
fi

############### add ROOTs for cmake of AxiSEM3D ###############
export MY_AXISEM_ROOTS=$HOME/.axisem3d_roots
if [ -f $MY_AXISEM_ROOTS ]; then
    rm $MY_AXISEM_ROOTS
fi
echo "export BOOST_ROOT=$MY_BOOST_DIR" >> $MY_AXISEM_ROOTS
echo "export FFTW_ROOT=$MY_FFTW_DIR" >> $MY_AXISEM_ROOTS
echo "export METIS_ROOT=$MY_METIS_DIR" >> $MY_AXISEM_ROOTS
echo "export EIGEN3_ROOT=$MY_EIGEN3_DIR" >> $MY_AXISEM_ROOTS
echo "export HDF5_ROOT=$MY_HDF5_DIR" >> $MY_AXISEM_ROOTS
if ! grep -Fxq ". $MY_AXISEM_ROOTS" $MY_BASH_PROFILE 
then
    echo ". $MY_AXISEM_ROOTS" >> $MY_BASH_PROFILE
fi
. $MY_BASH_PROFILE
