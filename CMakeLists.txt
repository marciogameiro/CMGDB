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
  # Python extension modules must use the dynamic CRT (/MD, CMake's MSVC
  # default): CPython itself is /MD, and mixing C runtimes across the
  # extension boundary is unsafe. Boost still links statically into the
  # module -- via its dynamic-runtime static libraries, so
  # Boost_USE_STATIC_RUNTIME must stay off to select them.
  set(Boost_USE_STATIC_LIBS ON)
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

# Vendored header-only sdsl-lite v3 (see src/CMGDB/_cmgdb/third_party/sdsl-lite/),
# compiled only when CMGDB_USE_SUCCINCT is defined (it backs the optional
# succinct grid; default builds use PointerGrid and never include sdsl).
# BEFORE, so that an sdsl v2 installation under a system prefix (e.g.
# /usr/local or Homebrew) cannot shadow these headers.
target_include_directories(_cmgdb BEFORE PRIVATE
  ${CMAKE_SOURCE_DIR}/${USER_PROJECT_PATH}/third_party/sdsl-lite/include)

target_link_libraries(_cmgdb PRIVATE Boost::chrono Boost::thread Boost::serialization)

install(TARGETS _cmgdb DESTINATION ${SKBUILD_PROJECT_NAME})
