# axisem3d_depends.sh
# created by Kuangdai on 30-Oct-2016 
# one-click installation of AxiSEM3D dependencies


############################## edit within this box ##############################
# installa path of brew
# Mac OS X default
export MY_BREW_INSTALL_DIR=/usr/local/
# Linux default
# export MY_BREW_INSTALL_DIR=$HOME/.linuxbrew/

# your bash_profile
# Mac OS X default
export MY_BASH_PROFILE=$HOME/.bash_profile
# Linux default
# export MY_BASH_PROFILE=$HOME/.bashrc

# What do you have already? And where are they?
# mpi
export MY_MPI_READY=false
# boost (make sure it is built with boost::mpi)
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
# NetCDF
export MY_NETCDF_READY=false
export MY_NETCDF_DIR=""
# exodus
export MY_EXODUS_READY=false
export MY_EXODUS_DIR=""

# specify where you want to install eigen3 and exodus
export EIGEN3_INSTALL_DIR=$HOME/axisem3d_depends/eigen3 
export EXODUS_INSTALL_DIR=$HOME/axisem3d_depends/exodus 
############################## edit within this box ##############################


##### cmake #####
brew install cmake

##### MPI #####
if [ $MY_MPI_READY == false ]; then
    brew install open-mpi
fi

##### Boost #####
if [ $MY_BOOST_READY == false ]; then
    brew install boost --with-mpi
    export MY_BOOST_DIR=$MY_BREW_INSTALL_DIR
    export MY_BOOST_READY=true
fi

##### FFTW #####
if [ $MY_FFTW_READY == false ]; then
    brew install fftw
    export MY_FFTW_DIR=$MY_BREW_INSTALL_DIR
    export MY_FFTW_READY=true
fi

##### METIS #####
if [ $MY_METIS_READY == false ]; then
    brew install homebrew/science/metis
    export MY_METIS_DIR=$MY_BREW_INSTALL_DIR
    export MY_METIS_READY=true
fi

##### EIGEN3 #####
if [ $MY_EIGEN3_READY == false ]; then
    # Brew does not have the correct version we need
    mkdir -p $EIGEN3_INSTALL_DIR
    if [ ! -f $EIGEN3_INSTALL_DIR/3.3-rc1.tar.bz2 ]; then
        wget -P $EIGEN3_INSTALL_DIR http://bitbucket.org/eigen/eigen/get/3.3-rc1.tar.bz2
    fi
    tar -jxvf $EIGEN3_INSTALL_DIR/3.3-rc1.tar.bz2 --strip 1 -C $EIGEN3_INSTALL_DIR
    export MY_EIGEN3_READY=true
    export MY_EIGEN3_DIR=$EIGEN3_INSTALL_DIR
fi

##### HDF5 #####
if [ $MY_HDF5_READY == false ]; then
    brew install homebrew/science/hdf5
    export MY_HDF5_DIR=$MY_BREW_INSTALL_DIR
    export MY_HDF5_READY=true
fi

##### NetCDF #####
if [ $MY_NETCDF_READY == false ]; then
    brew install homebrew/science/netcdf
    export MY_NETCDF_DIR=$MY_BREW_INSTALL_DIR
    export MY_NETCDF_READY=true
fi

##### Exodus #####
if [ $MY_EXODUS_READY == false ]; then
    if [ ! -f $EXODUS_INSTALL_DIR/exodus.zip ]; then
        wget -P $EXODUS_INSTALL_DIR https://github.com/gsjaardema/seacas/archive/exodus.zip
    fi    
    unzip -o $EXODUS_INSTALL_DIR/exodus.zip -d $EXODUS_INSTALL_DIR
    # tell cmake-exodus hdf5 and netcdf paths
    perl -pi -w -e '!$x && s/{TPL}/{MY_NETCDF_DIR}/ && ($x=1)' $EXODUS_INSTALL_DIR/seacas-exodus/cmake-exodus
    perl -pi -w -e 's/{TPL}/{MY_HDF5_DIR}/' $EXODUS_INSTALL_DIR/seacas-exodus/cmake-exodus
    # build
    my_start_point=$PWD
    mkdir $EXODUS_INSTALL_DIR/seacas-exodus/build
    cd $EXODUS_INSTALL_DIR/seacas-exodus/build
    ../cmake-exodus
    make -j4
    make install
    cd $my_start_point
    export MY_EXODUS_READY=true
    export MY_EXODUS_DIR=$EXODUS_INSTALL_DIR/seacas-exodus
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
echo "export EXODUS_ROOT=$MY_EXODUS_DIR" >> $MY_AXISEM_ROOTS
if ! grep -Fxq ". $MY_AXISEM_ROOTS" $MY_BASH_PROFILE 
then
    echo ". $MY_AXISEM_ROOTS" >> $MY_BASH_PROFILE
fi
. $MY_BASH_PROFILE
