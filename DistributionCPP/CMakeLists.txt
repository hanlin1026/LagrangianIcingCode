cmake_minimum_required(VERSION 2.6)
project(icing)

add_definitions( -std=c++0x -c -g )

include_directories( ${CMAKE_SOURCE_DIR} )
add_library( IcingLib
  Grid/PLOT3D.cc
  QuadTree/Bucket.cpp
  Cloud/Cloud.cpp
  Cloud/State.cpp )


add_executable( IcingDriver IcingDriver.cpp )
target_link_libraries( IcingDriver IcingLib )