cmake_minimum_required(VERSION 2.6)
project(icing)

add_definitions( -std=c++0x '-DCOMPLEX=std::complex<double>' -c -g -lmv -lsparse -lspblas )

include_directories( ${CMAKE_SOURCE_DIR}
		     /usr/include/eigen3/
		     /usr/include/gsl/
		     /usr/lib/SparseLib++/1.7/include/
                     /usr/lib/SparseLib++/1.7/mv/include/
 )
add_library( IcingLib
  Grid/PLOT3D.cc
  QuadTree/Bucket.cpp
  Cloud/Cloud.cpp
  Cloud/State.cpp
  Cloud/calcImpingementLimits.cpp
  Airfoil/Airfoil.cpp
  InputData/readInputParams.cpp
  AutoGridGen/autoGridGen.cpp
  ThermoEqns/ThermoEqns.cpp
  findAll.cpp )


add_executable( CATFISH IcingDriver.cpp )
target_link_libraries( CATFISH 
                       IcingLib
                       /usr/lib/libgsl.a 
		       /usr/lib/SparseLib++/1.7/lib/libmv.a
		       /usr/lib/SparseLib++/1.7/lib/libsparse.a
		       /usr/lib/SparseLib++/1.7/lib/libspblas.a )