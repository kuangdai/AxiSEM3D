# CMakeLists.txt
# created by Kuangdai on 26-Mar-2016
# CMakeLists.txt for AxiSEM3D

cmake_minimum_required(VERSION 3.0.0)
project(AxiSEM3D)
enable_language(Fortran)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

##################################### EDITABLE #####################################
# Set all necessary variables here.
# Please only try to edit within this block.
# If something is missing add it.

# compiler suit
SET(CMAKE_C_COMPILER           mpicc)
SET(CMAKE_CXX_COMPILER         mpicxx)
SET(CMAKE_Fortran_COMPILER     gfortran)

# external packages
# uncomment and edit the following lines to specify preferred prefix
# SET(EIGEN3_ROOT           edit_here)
# SET(FFTW_ROOT             edit_here)
# SET(BOOST_ROOT            edit_here)
# SET(METIS_ROOT            edit_here)
# SET(NETCDF_ROOT           edit_here)
# SET(HDF5_ROOT             edit_here)

# polynomial order of spectral elements (from 1 to 8)
SET(NPOL 4)

# solver precision
SET(USE_DOUBLE FALSE)

# use parallel NetCDF or not
# Modules that require parallel NetCDF: 
# * seismic inversion
# * wavefield and model visualisation
SET(USE_PARALLEL_NETCDF FALSE)

# add dynamic link to HDF5 libraries, depending on your netcdf build
# if you have netcdf installed from anaconda, you won't need this usually
# if linking of axisem3d fails, set this to TRUE and specify HDF5_ROOT
# NOTE: HDF5_ROOT should point to the hdf5 library you used to build netcdf
SET(NEED_HDF5_LINK FALSE)

# additional libraries to link with
# SET(ADDITIONAL_LIBS "-lcurl")

# directory to store FFTW wisdom files
# Specify any directory that does not require a `sudo` to write.
# Once set, it is unlikely to be changed. 
SET(FFTW_WISDOM_DIR $ENV{HOME}/.axisem3d_fftw_wisdom_dir)
##################################### EDITABLE #####################################



############# macros used in solver #############
# NPOL
ADD_DEFINITIONS(-D_NPOL=${NPOL})

# USE_DOUBLE
if (USE_DOUBLE)
    ADD_DEFINITIONS(-D_USE_DOUBLE)
endif ()

# FFTW wisdom dir
ADD_DEFINITIONS(-D_FFTW_WISDOM_DIR=\"${FFTW_WISDOM_DIR}\")

# SOURCE DIR
ADD_DEFINITIONS(-D_PROJECT_DIR=\"${PROJECT_SOURCE_DIR}\")

# measure all costs in time loop
ADD_DEFINITIONS(-D_MEASURE_TIMELOOP)

# serial build
SET(SERIAL_BUILD FALSE)
if (SERIAL_BUILD)
    ADD_DEFINITIONS(-D_SERIAL_BUILD)
    SET(USE_PARALLEL_NETCDF FALSE)
endif ()

# parallel NetCDF
if (USE_PARALLEL_NETCDF)
    ADD_DEFINITIONS(-D_USE_PARALLEL_NETCDF)
endif ()

############# find packages #############
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# mpi
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
# eigen
find_package(EIGEN3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
# fftw
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIR})
# boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
# metis
find_package(METIS REQUIRED)
include_directories(${METIS_INCLUDE_DIR})
# netcdf
find_package(NETCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIR})
# hdf5
if (NEED_HDF5_LINK)
    if (HDF5_ROOT)
        set(ENV{HDF5_ROOT} ${HDF5_ROOT})
    endif ()
    find_package(HDF5 COMPONENTS C COMPONENTS HL REQUIRED)
endif ()

############# local includes #############
include_directories(
    src
    ############################## core ##############################
    src/core
    src/core/fftw

    src/core/point
    src/core/point/mass
    src/core/point/solid_fluid

    src/core/element
    src/core/element/crd
    src/core/element/fft
    src/core/element/grad
    src/core/element/prt

    src/core/element/material/attenuation
    src/core/element/material/attenuation/1D
    src/core/element/material/attenuation/3D
    src/core/element/material/elastic
    src/core/element/material/elastic/1D
    src/core/element/material/elastic/3D
    src/core/element/material/acoustic

    src/core/source
    src/core/output/pointwise
    src/core/output/surface
    src/core/domain
    src/core/newmark

    ############################## preloop ##############################
    src/preloop
    src/preloop/utilities
    src/preloop/utilities/netcdf
    src/preloop/spectral
    src/preloop/spectral/mapping

    src/preloop/exodus
    src/preloop/graph
    src/preloop/mesh

    src/preloop/nrfield
    src/preloop/nrfield/wisdom

    src/preloop/physics/attenuation
    src/preloop/physics/material
    src/preloop/physics/relabelling

    src/preloop/source
    src/preloop/source/stf
    src/preloop/source/offaxis
    src/preloop/receiver

    src/3d_model/3d_volumetric
    src/3d_model/3d_volumetric/s20_s40rts
    src/3d_model/3d_volumetric/crust1
    src/3d_model/3d_volumetric/simple_shapes
    src/3d_model/3d_volumetric/EMC
    src/3d_model/3d_geometric
    src/3d_model/3d_geometric/ellipticity
    src/3d_model/3d_geometric/crust1
    src/3d_model/3d_geometric/EMC
    src/3d_model/3d_oceanload
    src/3d_model/3d_oceanload/const
    src/3d_model/3d_oceanload/crust1
)

############# local source #############
add_executable(
    axisem3d
    src/main.cpp
    src/axisem.cpp
    src/ftz.c

    ############################## core ##############################
    src/core/fftw/SolverFFTW.cpp
    src/core/fftw/SolverFFTW_1.cpp
    src/core/fftw/SolverFFTW_3.cpp
    src/core/fftw/SolverFFTW_N3.cpp
    src/core/fftw/SolverFFTW_N6.cpp
    src/core/fftw/SolverFFTW_N9.cpp

    src/core/point/mass/Mass1D.cpp
    src/core/point/mass/Mass3D.cpp
    src/core/point/mass/MassOcean1D.cpp
    src/core/point/mass/MassOcean3D.cpp

    src/core/point/Point.cpp
    src/core/point/FluidPoint.cpp
    src/core/point/SolidPoint.cpp
    src/core/point/SolidFluidPoint.cpp
    src/core/point/solid_fluid/SFCoupling1D.cpp
    src/core/point/solid_fluid/SFCoupling3D.cpp

    src/core/element/material/acoustic/Acoustic1D.cpp
    src/core/element/material/acoustic/Acoustic3D.cpp

    src/core/element/material/attenuation/1D/Attenuation1D_Full.cpp
    src/core/element/material/attenuation/1D/Attenuation1D_CG4.cpp
    src/core/element/material/attenuation/3D/Attenuation3D_Full.cpp
    src/core/element/material/attenuation/3D/Attenuation3D_CG4.cpp

    src/core/element/material/elastic/1D/Elastic1D.cpp
    src/core/element/material/elastic/1D/Isotropic1D.cpp
    src/core/element/material/elastic/1D/TransverselyIsotropic1D.cpp
    src/core/element/material/elastic/1D/Anisotropic1D.cpp
    src/core/element/material/elastic/3D/Elastic3D.cpp
    src/core/element/material/elastic/3D/Isotropic3D.cpp
    src/core/element/material/elastic/3D/TransverselyIsotropic3D.cpp
    src/core/element/material/elastic/3D/Anisotropic3D.cpp

    src/core/element/grad/Gradient.cpp
    
    src/core/element/crd/CrdTransTIsoSolid.cpp
    src/core/element/crd/CrdTransTIsoFluid.cpp
    src/core/element/fft/FieldFFT.cpp
    src/core/element/prt/PRT_1D.cpp
    src/core/element/prt/PRT_3D.cpp
    
    src/core/element/Element.cpp
    src/core/element/SolidElement.cpp
    src/core/element/FluidElement.cpp

    src/core/source/SourceTerm.cpp
    src/core/source/SourceTimeFunction.cpp
    src/core/output/pointwise/PointwiseRecorder.cpp
    src/core/output/pointwise/PointwiseIOAscii.cpp
    src/core/output/pointwise/PointwiseIONetCDF.cpp
    src/core/output/surface/SurfaceRecorder.cpp
    src/core/output/surface/SurfaceIO.cpp
    src/core/output/surface/SurfaceInfo.cpp
    src/core/domain/Domain.cpp
    src/core/newmark/Newmark.cpp

    ############################## preloop ##############################
    src/preloop/utilities/XMath.cpp
    src/preloop/utilities/Geodesy.cpp
    src/preloop/utilities/XMPI.cpp
    src/preloop/utilities/Parameters.cpp
    src/preloop/utilities/PreloopGradient.cpp
    src/preloop/utilities/PreloopFFTW.cpp
    src/preloop/utilities/MultilevelTimer.cpp
    src/preloop/utilities/netcdf/NetCDF_Reader.cpp
    src/preloop/utilities/netcdf/NetCDF_ReaderAscii.cpp
    src/preloop/utilities/netcdf/NetCDF_Writer.cpp

    src/preloop/spectral/SpectralConstants.cpp
    src/preloop/spectral/mapping/Mapping.cpp
    src/preloop/spectral/mapping/SphericalMapping.cpp
    src/preloop/spectral/mapping/LinearMapping.cpp
    src/preloop/spectral/mapping/SemiSphericalMapping.cpp

    src/preloop/exodus/ExodusModel.cpp
    src/preloop/graph/DualGraph.cpp
    src/preloop/graph/Connectivity.cpp
    src/preloop/mesh/GLLPoint.cpp
    src/preloop/mesh/Quad.cpp
    src/preloop/mesh/Mesh.cpp
    src/preloop/mesh/SlicePlot.cpp

    src/preloop/nrfield/NrField.cpp
    src/preloop/nrfield/ConstNrField.cpp
    src/preloop/nrfield/EmpNrField.cpp
    src/preloop/nrfield/UserNrField.cpp
    src/preloop/nrfield/wisdom/WisdomNrField.cpp
    src/preloop/nrfield/wisdom/NuWisdom.cpp

    src/preloop/physics/attenuation/AttBuilder.cpp
    src/preloop/physics/attenuation/AttSimplex.cpp
    src/preloop/physics/attenuation/AttAxiSEM.cpp
    src/preloop/physics/attenuation/simplex.f90

    src/preloop/physics/material/Material.cpp
    src/preloop/physics/relabelling/Relabelling.cpp

    src/preloop/source/Source.cpp
    src/preloop/source/Earthquake.cpp
    src/preloop/source/PointForce.cpp
    src/preloop/source/NullSource.cpp
    src/preloop/source/stf/STF.cpp
    src/preloop/source/stf/ErfSTF.cpp
    src/preloop/source/stf/GaussSTF.cpp
    src/preloop/source/stf/RickerSTF.cpp
    src/preloop/source/offaxis/OffAxisSource.cpp
    src/preloop/source/offaxis/OffAxisPointForce.cpp
    src/preloop/receiver/Receiver.cpp
    src/preloop/receiver/ReceiverCollection.cpp

    src/3d_model/3d_volumetric/Volumetric3D.cpp
    src/3d_model/3d_volumetric/s20_s40rts/s20rts.f90
    src/3d_model/3d_volumetric/s20_s40rts/s40rts.f90
    src/3d_model/3d_volumetric/s20_s40rts/Volumetric3D_s20rts.cpp
    src/3d_model/3d_volumetric/s20_s40rts/Volumetric3D_s40rts.cpp
    src/3d_model/3d_volumetric/crust1/Volumetric3D_crust1.cpp
    src/3d_model/3d_volumetric/simple_shapes/Volumetric3D_bubble.cpp
    src/3d_model/3d_volumetric/simple_shapes/Volumetric3D_cylinder.cpp
    src/3d_model/3d_volumetric/EMC/Volumetric3D_EMC.cpp

    src/3d_model/3d_geometric/Geometric3D.cpp
    src/3d_model/3d_geometric/ellipticity/Ellipticity.cpp
    src/3d_model/3d_geometric/crust1/Geometric3D_crust1.cpp
    src/3d_model/3d_geometric/EMC/Geometric3D_EMC.cpp

    src/3d_model/3d_oceanload/OceanLoad3D.cpp
    src/3d_model/3d_oceanload/const/OceanLoad3D_const.cpp
    src/3d_model/3d_oceanload/crust1/OceanLoad3D_crust1.cpp
)

############# link #############
target_link_libraries(
    axisem3d
    ${MPI_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${METIS_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${ADDITIONAL_LIBS}
)
