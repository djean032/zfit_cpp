cmake_minimum_required(VERSION 3.5)

project(zfit_cpp)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O3")

find_package(Ceres REQUIRED)
find_package(SUNDIALS REQUIRED)

set(EXT_INSTALL_PATH /usr/local)

set(SUNDIALS_INCLUDE_DIR
  "${EXT_INSTALL_PATH}/include"
  CACHE PATH "Location of SUNDIALS header files")

if(EXISTS "${EXT_INSTALL_PATH}/lib")
  set(SUNDIALS_LIBRARY_DIR
    "${EXT_INSTALL_PATH}/lib"
    CACHE PATH "Location of SUNDIALS libraries")
else()
  set(SUNDIALS_LIBRARY_DIR
    "${EXT_INSTALL_PATH}/lib64"
    CACHE PATH "Location of SUNDIALS libraries")
endif()
message(STATUS "SUNDIALS Library Dir: ${SUNDIALS_LIBRARY_DIR}")

# Find the SUNDIALS libraries
find_library(SUNDIALS_SOLVER_LIB
  sundials_cvode ${SUNDIALS_LIBRARY_DIR}
  DOC "CVODE library")

find_library(SUNDIALS_CORE_LIB
  sundials_core ${SUNDIALS_LIBRARY_DIR}
  DOC "SUNDIALS core library")

find_library(SUNDIALS_NVEC_LIB
  sundials_nvecserial ${SUNDIALS_LIBRARY_DIR}
  DOC "NVECTOR_SERIAL library")

find_library(SUNDIALS_MANYVEC_LIB
  sundials_nvecmanyvector ${SUNDIALS_LIBRARY_DIR}
  DOC "NVECTOR_MANYVECTOR library")

set(SUNDIALS_EXTRA_LIBS -lm CACHE STRING "Extra libraries needed by SUNDIALS")
#print out include directories
message(STATUS "SUNDIALS Include Dir: ${SUNDIALS_INCLUDE_DIR}")

set(SUNDIALS_LIBRARIES
  -L${SUNDIALS_LIBRARY_DIR}
  ${SUNDIALS_CORE_LIB}
  ${SUNDIALS_SOLVER_LIB}
  ${SUNDIALS_NVEC_LIB}
  ${SUNDIALS_MANYVEC_LIB}
  ${SUNDIALS_EXTRA_LIBS})

message(STATUS "SUNDIALS Libraries: ${SUNDIALS_LIBRARIES}")

find_program(CLANG_FORMAT_EXECUTABLE clang-format)
if(CLANG_FORMAT_EXECUTABLE)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXECUTABLE} -i ${CMAKE_SOURCE_DIR}
    COMMENT "Running clang-format"
  )
endif()

# helloworld
file(GLOB_RECURSE EDV_SOURCE_FILES "src/*.cpp" "src/*.hpp")
add_library(${PROJECT_NAME} ${EDV_SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC ${SUNDIALS_LIBRARIES})
add_executable(zfit zfit.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${SUNDIALS_INCLUDE_DIR})
target_link_libraries(zfit ${PROJECT_NAME})
