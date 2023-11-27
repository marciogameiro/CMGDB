cmake_minimum_required(VERSION 3.13)
# cmake_minimum_required(VERSION 2.8.12)

project(CMGDB)

cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0127 NEW)

# Visual Studio C++ does not support keywords such
# as "and", "not", etc. Settting the /permissive-
# flag for the compiler makes it support them.
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /permissive-")
endif()

if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
endif(APPLE)

add_subdirectory(pybind11)

find_package(Boost COMPONENTS chrono thread serialization REQUIRED)

include_directories (
  /usr/local/include
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/include
  ${USER_INCLUDE_PATH}/database
  ${USER_INCLUDE_PATH} )

link_directories (
  /usr/local/lib
  ${USER_LIBRARY_PATH} )

pybind11_add_module(_cmgdb src/CMGDB/_cmgdb/CMGDB.cpp)

target_link_libraries(_cmgdb PRIVATE Boost::chrono Boost::thread Boost::serialization sdsl)
