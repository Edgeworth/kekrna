cmake_minimum_required(VERSION 3.2)

project(kekrna)

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
file(GLOB KEKRNA_SOURCE "src/*.cpp")
file(GLOB KEKRNA_HEADERS "src/*.h")
list(REMOVE_ITEM KEKRNA_SOURCE src/main.cpp)
file(GLOB TEST_SOURCE "tests/*.cpp")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindUBSan)
include(FindASan)
include(FindMSan)
include(FindThreads)

include_directories(src)

add_library(kekrna_lib ${KEKRNA_SOURCE} ${KEKRNA_HEADERS})
add_executable(kekrna src/main.cpp)
target_link_libraries(kekrna kekrna_lib)

add_subdirectory(lib/gtest)
enable_testing()
add_executable(run_tests ${TEST_SOURCE} ${KEKRNA_HEADERS})
target_link_libraries(run_tests Threads::Threads gtest kekrna_lib)
