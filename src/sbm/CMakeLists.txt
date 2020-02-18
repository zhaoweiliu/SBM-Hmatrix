cmake_minimum_required(VERSION 3.2)
project(sbm C CXX)

set(SOURCE_FILES
    division.h
    division.cpp
)

add_library(sbm SHARED STATIC ${SOURCE_FILES})

install(TARGETS sbm DESTINATION ${SBM_INSTALL_LIB_DIR})
install(FILES division.h DESTINATION ${SBM_INSTALL_INCLUDE_DIR})