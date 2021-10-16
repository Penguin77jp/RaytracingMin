cmake_minimum_required(VERSION 3.20)
project(RaytracingMin)

set(CMAKE_CXX_STANDARD 20)

include_directories(RaytracingMin json/include stb src/*.h)

file(GLOB CXX_FILES src/*.cpp)
add_executable(RaytracingMin src/main.cpp ${CXX_FILES})
