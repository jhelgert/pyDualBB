cmake_minimum_required(VERSION 3.16.0)
project(dualBB VERSION 0.1.0)
include(FetchContent)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx2 -DNDEBUG")

set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx2 -DNDEBUG -fPIC")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Cplex REQUIRED)
include(Blaze)
include(nlohmann-json)
include(pybind11)

add_library(Helpers STATIC ${SOURCE_DIR}/Helpers.cpp ${SOURCE_DIR}/Helpers.hpp)
target_include_directories(Helpers PUBLIC ${blaze_SOURCE_DIR})

add_library(cplexDualBB STATIC ${SOURCE_DIR}/cplexDualBB.cpp ${SOURCE_DIR}/cplexDualBB.hpp)
target_include_directories(cplexDualBB PUBLIC ${CPLEX_INCLUDE_DIRS})
target_include_directories(cplexDualBB PUBLIC ${blaze_SOURCE_DIR})
target_link_libraries(cplexDualBB PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(cplexDualBB PUBLIC ${BLAS_LIBRARIES})
target_link_libraries(cplexDualBB PUBLIC ${LAPACK_LIBRARIES})
target_link_libraries(cplexDualBB PUBLIC ${CPLEX_LIBRARIES})

add_library(cplexLogging STATIC ${SOURCE_DIR}/cplexLogging.cpp ${SOURCE_DIR}/cplexLogging.hpp)
target_include_directories(cplexLogging PUBLIC ${CPLEX_INCLUDE_DIRS})
target_link_libraries(cplexLogging PUBLIC ${CPLEX_LIBRARIES})

pybind11_add_module(_cplex_dualbb_wrapper ${SOURCE_DIR}/python_binding.cpp)
target_include_directories(_cplex_dualbb_wrapper PUBLIC ${CPLEX_INCLUDE_DIRS})
target_include_directories(_cplex_dualbb_wrapper PUBLIC ${blaze_SOURCE_DIR})
target_link_libraries(_cplex_dualbb_wrapper PUBLIC Helpers cplexDualBB cplexLogging)
target_link_libraries(_cplex_dualbb_wrapper PUBLIC ${BLAS_LIBRARIES})
target_link_libraries(_cplex_dualbb_wrapper PUBLIC ${LAPACK_LIBRARIES})
target_link_libraries(_cplex_dualbb_wrapper PRIVATE nlohmann_json::nlohmann_json)

