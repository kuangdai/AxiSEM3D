# AxiSEM3D - A Quickstart Guide

I look ugly to you? I am a Markdown file. If you do not have a Markdown reader, read me on [github](https://github.com/kuangdai/AxiSEM3D), or paste me [here](http://dillinger.io/).

Report issues to kuangdal@earth.ox.ac.uk.

### Show-off videos on Youtube:
[![](https://img.youtube.com/vi/60cMKohNQoU/2.jpg)](https://www.youtube.com/watch?v=60cMKohNQoU)
[![](https://img.youtube.com/vi/nw98Xxy4TdM/2.jpg)](https://www.youtube.com/watch?v=nw98Xxy4TdM)

## 1 Get AxiSEM3D
```sh
git clone https://github.com/kuangdai/AxiSEM3D
export CURRENT_WORK_DIR=$PWD
```

## 2 Installing dependencies
AxiSEM3D has been built upon a few modern numerical packages for its performance and sustainability, as listed below (newer versions are acceptable):

name | version | build-from-source instructions
--- | --- | --- 
MPI | --- | Common implementations such as [open-mpi](https://www.open-mpi.org/) or [mpich2](http://www.mpich.org/).
Boost | 1.60 | Simply [download](https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.bz2) and unzip to your path of treasures. No need to build or install.
Eigen3 | 3.3-rc1 | Simply [download](http://bitbucket.org/eigen/eigen/get/3.3-rc1.tar.bz2) and unzip to your path of treasures. No need to build or install.
FFTW | 3.3 | Install both double-precision (without _--enable-float_) and single-precision (with _--enable-float_) versions to the same path. See [here](http://www.fftw.org/fftw2_doc/fftw_6.html).
METIS | 5.1 | [Download](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) and install following `Install.txt`. Skip step 3 in `Install.txt` to use the 32-bit build. 
NetCDF | 4.4 | Follow instructions [here](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html). To avoid compatibility issues, AxiSEM3D does not directly use HDF5, the major dependency of NetCDF4. 

 Don't panic! All these popular packages may be handily installed with free package management software, such as [Conda](http://conda.pydata.org/docs/). Here we introduce the wizard `axisem3d_depends.sh`. 

* Make sure your MPI works properly. 
* Get [Conda](http://conda.pydata.org/docs/).  
* Edit the first few lines in `axisem3d_depends.sh` and run it.
* Check your `~/.bash_profile` (or `~/.bashrc`) and `~/.axisem3d_roots`. Grats, you are done!
    

## 3 Building AxiSEM3D
* Edit `SOLVER/CMakeLists.txt` if needed (normally not), including

    to-be-edited | notes
    ---|---
    compiler suit | Changes are usually required only for HPC clusters, unless you have a bizarre MPI. 
    dependency roots | No need if you have used `axisem3d_depends.sh` to install the dependencies.
    FFTW_WISDOM_DIR | Just specify any directory you like, or leave it as it. 
 
* Build AxiSEM3D (go step by step to see what's happening)

    * Style 1: work under source 

        ```sh
        cd $CURRENT_WORK_DIR/AxiSEM3D
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
      -- Build files have been written to: ...
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
        cd $CURRENT_WORK_DIR
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

## 4 The MESHER
In the above examples, we use the mesh file `AxiSEM_prem_ani_one_crust_50.e` (anisotropic PREM model with one crustal layer at a 50 s period), located at `SOLVER/template/input`. 

To generate an AxiSEM3D mesh like this, you will need the `salvus_mesher_lite`, a python-based command-line tool to generate AxiSEM/AxiSEM3D meshes (credit to van Driel Martin, Lion Krischer and others from the Salvus group at ETH Zurich). 

* See installation and usage of `salvus_mesher_lite` from [here](https://gitlab.com/Salvus/SalvusMesherLite). 



