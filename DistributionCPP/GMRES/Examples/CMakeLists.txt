cmake_minimum_required(VERSION 2.6)
project(gmres)

add_definitions('-DCOMPLEX=std::complex<double>' -c -g -lmv -lsparse -lspblas)

include_directories( ${CMAKE_SOURCE_DIR}
                     /usr/lib/SparseLib++/1.7/include/
                     /usr/lib/SparseLib++/1.7/mv/include/)

add_executable( GMRES tgmres.cc )

target_link_libraries( GMRES
                       /usr/lib/SparseLib++/1.7/lib/libmv.a
		       /usr/lib/SparseLib++/1.7/lib/libsparse.a
		       /usr/lib/SparseLib++/1.7/lib/libspblas.a)