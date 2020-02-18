cmake_minimum_required(VERSION 3.2)
project(sbm_tests)

add_subdirectory(lib/googletest)

include_directories(${SBM_HEADERS_DIR})
include_directories(lib/googletest/googletest/include)

set(SOURCE_FILES main.cpp src/divider_tests.cpp)

add_executable(sbm_tests ${SOURCE_FILES})

target_link_libraries(sbm_tests sbm gtest)

install(TARGETS sbm_tests DESTINATION bin)

