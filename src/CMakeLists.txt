cmake_minimum_required(VERSION 3.2)
project(sbm_hmatrix)

add_subdirectory(sbm)
set(SOURCE_FILES main.cpp)

add_executable(sbm_hmatrix ${SOURCE_FILES})

target_link_libraries(sbm_hmatrix sbm hpro)

install(TARGETS sbm_hmatrix DESTINATION ${SBM_INSTALL_BIN_DIR})
