cmake_minimum_required(VERSION 3.13)

project(CMGDB)

cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0127 NEW)

# Visual Studio C++ does not support keywords such
# as "and", "not", etc. Settting the /permissive-
# flag for the compiler makes it support them.
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /permissive-")
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /permissive-")
  set(Boost_USE_STATIC_LIBS ON)
  set(Boost_USE_STATIC_RUNTIME ON)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  # Set compiler flag for sdsl-lite
  add_definitions("/DMSVC_COMPILER")
endif()

if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
endif(APPLE)

add_subdirectory(pybind11)

# set(CMAKE_TOOLCHAIN_FILE C:/vcpkg/scripts/buildsystems/vcpkg.cmake)

find_package(Boost COMPONENTS chrono thread serialization REQUIRED)

if(Boost_FOUND)
  message(STATUS "Boost version: ${Boost_VERSION}")
  message("BOOST_ROOT is ${BOOST_ROOT}")
  message("Boost_INCLUDE_DIRS is ${Boost_INCLUDE_DIRS}")
  message("Boost_LIBRARY_DIRS is ${Boost_LIBRARY_DIRS}")
  message("Boost_LIBRARIES is ${Boost_LIBRARIES}")
endif()

# Set Boost include and lib directories
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

if(MSVC)
  # Set GMP and SDSL include and link directories
  message("GMP_ROOT is $ENV{GMP_ROOT}")
  include_directories($ENV{GMP_ROOT}/include)
  link_directories ($ENV{GMP_ROOT}/lib)
  message("SDSL_ROOT is $ENV{SDSL_ROOT}")
  include_directories($ENV{SDSL_ROOT}/include)
  link_directories ($ENV{SDSL_ROOT}/lib/Release)
endif()

include_directories (
  ${USER_INCLUDE_PATH}/database
  ${USER_INCLUDE_PATH}
  /usr/local/include )

link_directories (
  ${USER_LIBRARY_PATH}
  /usr/local/lib )

pybind11_add_module(_cmgdb src/CMGDB/_cmgdb/CMGDB.cpp)

target_link_libraries(_cmgdb PRIVATE Boost::chrono Boost::thread Boost::serialization sdsl)
