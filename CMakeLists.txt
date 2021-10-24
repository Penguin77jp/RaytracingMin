cmake_minimum_required(VERSION 3.20)
project(RaytracingMin)

set(CMAKE_CXX_STANDARD 17)

# enable openmp
if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    add_compile_options(-D _DEBUG)
endif()
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(RaytracingMin json/include stb src/*.h)
file(GLOB CXX_FILES src/*.cpp)
add_executable(RaytracingMin ${CXX_FILES})
