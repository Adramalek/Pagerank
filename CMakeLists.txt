cmake_minimum_required(VERSION 3.5)
project(Pagerank)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(SOURCE_FILES
    src/main.cpp
    src/pagerank.cpp
    src/pagerank.h
    src/utils.cpp
    src/utils.h src/IPagerankMatrixReader.h)

add_executable(Pagerank ${SOURCE_FILES})