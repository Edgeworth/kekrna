cmake_minimum_required(VERSION 3.2)

project(kekrna)

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O2")
file(GLOB KEKRNA_SOURCE "src/*.cpp" "src/*.h")
file(GLOB TEST_SOURCE "tests/*.cpp" "tests/*.h")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindUBSan)
include(FindASan)
include(FindMSan)
include(FindThreads)

include_directories(src)

add_executable(compute_energy ${KEKRNA_SOURCE} src/programs/compute_energy.cpp)
add_executable(batch_energy ${KEKRNA_SOURCE} src/programs/batch_energy.cpp)
add_executable(fold ${KEKRNA_SOURCE} src/programs/fold.cpp)

add_subdirectory(lib/gtest)
enable_testing()
add_executable(run_tests ${TEST_SOURCE} ${KEKRNA_SOURCE})
target_link_libraries(run_tests Threads::Threads gtest)
