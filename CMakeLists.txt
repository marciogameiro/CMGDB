cmake_minimum_required(VERSION 3.15...3.26)

project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES CXX)

cmake_policy(SET CMP0144 NEW)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Note on optimization flags: the Release build already compiles with -O3.
# LTO was tried and measured to have no effect -- the project builds as a
# single translation unit (CMGDB.cpp includes the header-only codebase), so
# there are no cross-TU boundaries for LTO to optimize.

# Visual Studio C++ does not support keywords such
# as "and", "not", etc. Settting the /permissive-
# flag for the compiler makes it support them.
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /permissive-")
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /permissive-")
  set(Boost_USE_STATIC_LIBS ON)
  set(Boost_USE_STATIC_RUNTIME ON)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
endif()

if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 CONFIG REQUIRED)

# With the settings above pybind11 finds Python
# find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

find_package(Boost COMPONENTS chrono thread serialization CONFIG REQUIRED)

if(Boost_FOUND)
  message(STATUS "Boost version: ${Boost_VERSION}")
  message("Boost_ROOT is ${Boost_ROOT}")
  message("Boost_INCLUDE_DIRS is ${Boost_INCLUDE_DIRS}")
  message("Boost_LIBRARY_DIRS is ${Boost_LIBRARY_DIRS}")
  message("Boost_LIBRARIES is ${Boost_LIBRARIES}")
endif()

# Set Boost include and lib directories
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

if(MSVC)
  # Set GMP include and link directories
  message("GMP_ROOT is $ENV{GMP_ROOT}")
  include_directories($ENV{GMP_ROOT}/include)
  link_directories($ENV{GMP_ROOT}/lib)
endif()

if(APPLE)
  # Homebrew directories for Apple silicon
  # Homebrew include and library dirs (for GMP)
  include_directories(/opt/homebrew/include)
  link_directories(/opt/homebrew/lib)
endif()

# USER_PROJECT_PATH is set in pyproject.toml
include_directories(${CMAKE_SOURCE_DIR}/${USER_PROJECT_PATH}/include)
include_directories(${CMAKE_SOURCE_DIR}/${USER_PROJECT_PATH}/include/database)
# include_directories(${CMAKE_SOURCE_DIR}/src/CMGDB/_cmgdb/include)
# include_directories(${CMAKE_SOURCE_DIR}/src/CMGDB/_cmgdb/include/database)

pybind11_add_module(_cmgdb MODULE ${USER_PROJECT_PATH}/CMGDB.cpp)
# pybind11_add_module(_cmgdb MODULE src/CMGDB/_cmgdb/CMGDB.cpp)

# Vendored header-only sdsl-lite v3 (see src/CMGDB/_cmgdb/third_party/sdsl-lite/).
# BEFORE, so that an sdsl v2 installation under a system prefix (e.g.
# /usr/local or Homebrew) cannot shadow these headers.
target_include_directories(_cmgdb BEFORE PRIVATE
  ${CMAKE_SOURCE_DIR}/${USER_PROJECT_PATH}/third_party/sdsl-lite/include)

target_link_libraries(_cmgdb PRIVATE Boost::chrono Boost::thread Boost::serialization)

install(TARGETS _cmgdb DESTINATION ${SKBUILD_PROJECT_NAME})
