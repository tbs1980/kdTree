cmake_minimum_required(VERSION 2.8)

project(KDTree CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wfatal-errors")

add_executable(testKDTree testKDTree.cpp)
add_executable(NDExample NDExample.cpp)
