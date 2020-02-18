cmake_minimum_required(VERSION 3.2)

project(sbm-hmatrix)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

set(SBM_INSTALL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(SBM_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(SBM_INSTALL_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

set(SBM_HEADERS_DIR ${PROJECT_SOURCE_DIR}/src/sbm)

set(H_LIB_PRO_DIR ${PROJECT_SOURCE_DIR}/hlibpro-2.8/lib)
set(TBB_DIRS ${PROJECT_SOURCE_DIR}/tbb)
set(TBB_INCLUDE_DIRS ${TBB_DIRS}/include)
set(TBB_LIB_DIRS ${TBB_DIRS}/lib)

include_directories(${SBM_INSTALL_INCLUDE_DIR})
include_directories(${SBM_HEADERS_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories("${PROJECT_SOURCE_DIR}/hlibpro-2.8/include/")
include_directories(${TBB_INCLUDE_DIRS})

link_directories(${H_LIB_PRO_DIR})
link_directories(${TBB_LIB_DIRS})

add_subdirectory(src)
add_subdirectory(test)
