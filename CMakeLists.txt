cmake_minimum_required(VERSION 3.20)

project(zfit_cpp LANGUAGES C CXX)

add_library(zfit_core
  src/kinetics.cpp
  src/laser.cpp
  src/model.cpp
)

target_compile_features(zfit_core PUBLIC cxx_std_23)
target_include_directories(zfit_core PUBLIC include)

option(ZFIT_ENABLE_CVODES "Build the SUNDIALS CVODES wrapper" on)
option(ZFIT_ENABLE_FITTING "Build the cminpack fitting layer" off)
option(ZFIT_BUILD_BENCHMARKS "Build solver benchmark executables" off)

# Headless measurement workup, independent of the solver and optimizer.
add_library(zfit_workup
  src/workup/covariance.cpp
  src/workup/dark.cpp
  src/workup/matching.cpp
  src/workup/prepare.cpp
  src/workup/ratio.cpp
)
target_compile_features(zfit_workup PUBLIC cxx_std_23)
target_include_directories(zfit_workup PUBLIC include)

if(ZFIT_ENABLE_CVODES)
  find_package(SUNDIALS CONFIG REQUIRED COMPONENTS cvodes nvecserial sunmatrixdense sunlinsoldense)

  add_library(zfit_cvodes
    src/cvodes_integrator.cpp
  )

  target_compile_features(zfit_cvodes PUBLIC cxx_std_23)
  target_include_directories(zfit_cvodes PUBLIC include)
  target_link_libraries(zfit_cvodes
    PRIVATE
      zfit_core
      SUNDIALS::cvodes
      SUNDIALS::nvecserial
      SUNDIALS::sunmatrixdense
      SUNDIALS::sunlinsoldense
  )
endif()

if(ZFIT_ENABLE_FITTING)
  if(NOT ZFIT_ENABLE_CVODES)
    message(FATAL_ERROR "CVODES is required for zfit_fit")
  endif()
  include(FetchContent)

  set(CMINPACK_PRECISION d CACHE STRING "" FORCE)
  set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

  FetchContent_Declare(
    cminpack
    GIT_REPOSITORY https://github.com/devernay/cminpack.git
    GIT_TAG v1.3.14
  )
  FetchContent_MakeAvailable(cminpack)

  set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 5.0.1
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(eigen)

  add_library(zfit_fit src/fit.cpp src/whitening.cpp)
  target_compile_features(zfit_fit PUBLIC cxx_std_23)
  target_include_directories(zfit_fit PUBLIC include)
  target_link_libraries(zfit_fit PUBLIC zfit_cvodes PRIVATE cminpack Eigen3::Eigen)
endif()

if(ZFIT_BUILD_BENCHMARKS)
  if(NOT ZFIT_ENABLE_CVODES)
    message(FATAL_ERROR "CVODES is required for zfit_cvodes_benchmark")
  endif()

  add_executable(zfit_cvodes_benchmark
    benchmarks/cvodes_curve_benchmark.cpp
  )
  target_compile_features(zfit_cvodes_benchmark PRIVATE cxx_std_23)
  target_link_libraries(zfit_cvodes_benchmark PRIVATE zfit_cvodes)
endif()

include(CTest)

if(BUILD_TESTING)
  include(FetchContent)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(Catch2)

  add_executable(zfit_tests
    tests/kinetics_tests.cpp
    tests/laser_tests.cpp
    tests/model_tests.cpp
    tests/workup/covariance_tests.cpp
    tests/workup/dark_tests.cpp
    tests/workup/matching_tests.cpp
    tests/workup/prepare_tests.cpp
    tests/workup/ratio_tests.cpp
  )
  target_link_libraries(zfit_tests PRIVATE zfit_core zfit_workup Catch2::Catch2WithMain)

  if(ZFIT_ENABLE_CVODES)
    target_sources(zfit_tests PRIVATE
      tests/cvodes_integrator_tests.cpp
      tests/legacy_parity_tests.cpp
    )
    target_link_libraries(zfit_tests PRIVATE zfit_cvodes)
  endif()

  if(ZFIT_ENABLE_FITTING)
    target_sources(zfit_tests PRIVATE
      tests/fit_tests.cpp
      tests/whitening_tests.cpp
    )
    target_link_libraries(zfit_tests PRIVATE zfit_fit)
  endif()

  target_compile_definitions(zfit_tests PRIVATE
    ZFIT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
  )

  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
  include(Catch)
  catch_discover_tests(zfit_tests)
endif()

find_program(CLANG_FORMAT_EXE NAMES clang-format)
if(CLANG_FORMAT_EXE)
  file(GLOB_RECURSE ALL_SOURCE_FILES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/*.cpp"
  )
  add_custom_target(format
    COMMAND "${CLANG_FORMAT_EXE}" -i --style=file ${ALL_SOURCE_FILES}
    COMMENT "Formatting first-party C++ files"
    VERBATIM
  )
endif()
