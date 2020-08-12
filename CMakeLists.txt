cmake_minimum_required(VERSION 3.1)

project(CMGDB)

add_subdirectory(lib/pybind11)

find_package(Boost REQUIRED COMPONENTS chrono thread serialization)

include_directories (
  /usr/local/include
  /usr/include
  /opt/local/include
  ${HOME}/.local/include
  ${USER_INCLUDE_PATH}/database
  ${USER_INCLUDE_PATH} )

link_directories (
  /usr/local/lib
  /usr/lib
  ${HOME}/.local/lib
  ${USER_LIBRARY_PATH} )

pybind11_add_module(_cmgdb src/CMGDB/_cmgdb/CMGDB.cpp)

target_link_libraries (
  _cmgdb PRIVATE
  Boost::chrono
  Boost::thread
  Boost::serialization
  sdsl )
