# AxiSEM3D - A Quickstart Guide

I look ugly to you? I am a Markdown file. If you do not have a Markdown reader, read me on [github](https://github.com/kuangdai/AxiSEM3D), or paste me [here](dillinger.io/).

Report issues to kuangdal@earth.ox.ac.uk.

## 1 Installing dependencies
AxiSEM3D has been built upon a few modern numerical packages for its performance and sustainability, as listed below (newer versions are acceptable):

name | version | build-from-source instructions
--- | --- | --- 
mpi | --- | Common implementations such as [open-mpi](https://www.open-mpi.org/) or [mpich2](www.mpich.org/) are recommended.
boost | 1.60 | Install [Boost C++ Libraries](http://www.boost.org/) with _boost::mpi_, following [here](http://www.boost.org/doc/libs/1_62_0/doc/html/mpi/getting_started.html).
eigen3 | 3.3-rc1 | Simply [download](http://bitbucket.org/eigen/eigen/get/3.3-rc1.tar.bz2) and unzip to your path of treasures.
fftw | 3.3 | Install both double- and single-precision (with _--enable-float_) versions to the same path. See [here](http://www.fftw.org/fftw2_doc/fftw_6.html).
metis | 5.1 | [Download](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) and install following `Install.txt`. Skip step 3 in `Install.txt` to use the 32-bit build. 
Exodus |---| Must be built from source, upon both [HDF5](https://support.hdfgroup.org/HDF5/) and [NetCDF4](http://www.unidata.ucar.edu/blogs/news/entry/netcdf-4-4-1). Refer to `axisem3d_depends.sh`.

 Don't panic! Most of these popular packages (except Exodus) may be handily installed with free package management software, such as [Homebrew](http://brew.sh/) for Mac OS X and [Linuxbrew](http://linuxbrew.sh/) for common Linux distributions. Here we provide a Brew-based wizard for Mac OS X and Linux users.

* Install Brew
    * [Homebrew](http://brew.sh/) for Mac OS X
    * [Linuxbrew](http://linuxbrew.sh/) for Linux
* Edit the first few lines in `axisem3d_depends.sh` and run it. 
* Check `~/.bash_profile` (or `~/.bashrc`) and `~/.axisem3d_roots`. You are done!
    

## 2 Building AxiSEM3D
* Get the source

    ```sh
    git clone https://github.com/kuangdai/AxiSEM3D
    export CURRENT_WORK_DIR=$PWD
    ```
* Edit `$AxiSEM3D_SOURCE/SOLVER/CMakeLists.txt` if needed (normally not), including

    to-be-edited | notes
    ---|---
    compiler suit | Changes are usually required only for HPC clusters, unless you have a bizarre MPI. 
    dependency roots | No need if you have used `axisem3d_depends.sh` to install the dependencies.
    FFTW_WISDOM_DIR | Just specify any directory you like, or leave it as it. 
 
* Build AxiSEM3D

    * Style 1: work under source 

        ```sh
        cd AxiSEM3D
        # make a simulation directory
        mkdir my_first_run
        cd my_first_run
        # cmake (*** see the notes below ***)
        cmake -DCMAKE_BUILD_TYPE=release ../SOLVER
        # compile and link
        make -j4
        # copy the input folder
        cp -R ../template/input ./
        # run it with any number of processors
        mpirun -np 4 ./axisem3d
        # check the outputs
        ls output
        ```
        
        Read the output of `cmake` carefully, which should end up with (If not, verify dependency installations as prompted)
        
      ```sh
      -- Configuring done
      -- Generating done
      -- Build files have been written to: $AxiSEM3D_BUILD
      ```
    
    * Style 2: keep source clean (suggested)
    
        To offer maximal flexibility for various infrastructures, AxiSEM3D is organized such that the directories of source, build and simulations are fully independent of one another. For exmaple:
        
        ```sh
        ########## build ##########
        cd $CURRENT_WORK_DIR
        # make a build directory
        mkdir my_axisem3d_build
        cd my_axisem3d_build
        # cmake
        cmake -DCMAKE_BUILD_TYPE=release ../AxiSEM3D/SOLVER
        # compile and link
        make -j4
        
        ########## run ##########
        # make a simulation directory
        mkdir my_second_run
        cd my_second_run
        # copy the executable 
        cp ../my_axisem3d_build/axisem3d ./
        # copy the input folder
        cp -R ../AxiSEM3D/template/input ./
        # Though optional, it is always a good practice 
        # to backup the source code for reproducibility:
        cp -r ../AxiSEM3D/SOLVER ./
        # run it with any number of processors
        mpirun -np 4 ./axisem3d
        # check the outputs
        ls output
        ```
        




