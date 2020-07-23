cmake_minimum_required(VERSION 2.6)

project(CMGDB)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ftemplate-depth=2048")

find_package(Boost REQUIRED COMPONENTS chrono thread serialization)

include_directories (${CMAKE_SOURCE_DIR}/include/database ${CMAKE_SOURCE_DIR}/include)

include_directories (/usr/local/include /opt/local/include /usr/include ${HOME}/.local/include)

link_directories (/usr/local/lib /usr/lib ${HOME}/.local/lib)

add_executable(CMGDB src/CMGDB.cpp)

target_link_libraries(CMGDB Boost::chrono Boost::thread Boost::serialization sdsl)
