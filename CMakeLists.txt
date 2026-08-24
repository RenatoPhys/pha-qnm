cmake_minimum_required(VERSION 3.28)
project(pha_qnm VERSION 0.1.0 LANGUAGES C CXX)

option(PHA_QNM_BUILD_TESTS "Build the validation suite" ON)
option(PHA_QNM_ENABLE_HDF5 "Enable native HDF5 posterior/result I/O" OFF)
option(PHA_QNM_ENABLE_NETWORK "Enable libcurl download and OpenSSL checksums" OFF)
option(PHA_QNM_ENABLE_LAPACK "Enable the dense generalized eigensolver" OFF)
option(PHA_QNM_ENABLE_SANITIZERS "Enable address/undefined sanitizers where supported" OFF)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CompilerWarnings)
include(Sanitizers)

add_library(pha_qnm_core
  src/model/potentials.cpp
  src/background/horizon_series.cpp
  src/background/integrator.cpp
  src/spectral/chebyshev.cpp
  src/perturbations/decoupled.cpp
  src/hydro/linear_modes.cpp
  src/io/posterior.cpp
)
target_include_directories(pha_qnm_core PUBLIC include)
pha_qnm_set_warnings(pha_qnm_core)
if(PHA_QNM_ENABLE_SANITIZERS)
  pha_qnm_enable_sanitizers(pha_qnm_core)
endif()

if(PHA_QNM_ENABLE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
  target_sources(pha_qnm_core PRIVATE src/io/posterior_hdf5.cpp)
  target_link_libraries(pha_qnm_core PRIVATE HDF5::HDF5)
  target_compile_definitions(pha_qnm_core PUBLIC PHA_QNM_HAS_HDF5=1)
endif()

if(PHA_QNM_ENABLE_NETWORK)
  find_package(CURL REQUIRED)
  find_package(OpenSSL REQUIRED COMPONENTS Crypto)
  target_sources(pha_qnm_core PRIVATE src/io/download.cpp)
  target_link_libraries(pha_qnm_core PRIVATE CURL::libcurl OpenSSL::Crypto)
  target_compile_definitions(pha_qnm_core PUBLIC PHA_QNM_HAS_NETWORK=1)
endif()

if(PHA_QNM_ENABLE_LAPACK)
  find_package(LAPACK REQUIRED)
  target_sources(pha_qnm_core PRIVATE src/spectral/dense_gevp.cpp)
  target_link_libraries(pha_qnm_core PRIVATE LAPACK::LAPACK)
  target_compile_definitions(pha_qnm_core PUBLIC PHA_QNM_HAS_LAPACK=1)
endif()

add_executable(pha-qnm apps/pha_qnm_main.cpp)
target_link_libraries(pha-qnm PRIVATE pha_qnm_core)
pha_qnm_set_warnings(pha-qnm)

if(PHA_QNM_BUILD_TESTS)
  enable_testing()
  add_executable(pha_qnm_tests
    tests/test_main.cpp
    tests/test_potentials.cpp
    tests/test_horizon_series.cpp
    tests/test_background.cpp
    tests/test_chebyshev.cpp
    tests/test_operators.cpp
    tests/test_hydro.cpp
  )
  if(PHA_QNM_ENABLE_LAPACK)
    target_sources(pha_qnm_tests PRIVATE tests/test_dense_gevp.cpp)
  endif()
  target_link_libraries(pha_qnm_tests PRIVATE pha_qnm_core)
  target_include_directories(pha_qnm_tests PRIVATE tests)
  pha_qnm_set_warnings(pha_qnm_tests)
  if(PHA_QNM_ENABLE_SANITIZERS)
    pha_qnm_enable_sanitizers(pha_qnm_tests)
  endif()
  add_test(NAME pha_qnm_validation COMMAND pha_qnm_tests)
endif()

install(TARGETS pha-qnm RUNTIME DESTINATION bin)
