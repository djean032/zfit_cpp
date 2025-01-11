cmake_minimum_required(VERSION 3.5)

project(zfit_cpp)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER icpx)
set(CMAKE_CXX_FLAGS "-O2")

find_package(Ceres REQUIRED)
find_package(SUNDIALS 7.1.0 REQUIRED)
# helloworld
file(GLOB_RECURSE EDV_SOURCE_FILES "src/*.cpp" "src/*.hpp")
add_library(${PROJECT_NAME} ${EDV_SOURCE_FILES})
add_executable(zfit zfit.cpp)
target_link_libraries(zfit Ceres::ceres SUNDIALS::cvode ${PROJECT_NAME})
