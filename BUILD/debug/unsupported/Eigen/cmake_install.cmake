# Install script for directory: F:/Course/fem/eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/Eigen3")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "C:/mingw64/bin/objdump.exe")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "F:/Course/fem/eigen/unsupported/Eigen/AdolcForward"
    "F:/Course/fem/eigen/unsupported/Eigen/AlignedVector3"
    "F:/Course/fem/eigen/unsupported/Eigen/ArpackSupport"
    "F:/Course/fem/eigen/unsupported/Eigen/AutoDiff"
    "F:/Course/fem/eigen/unsupported/Eigen/BVH"
    "F:/Course/fem/eigen/unsupported/Eigen/EulerAngles"
    "F:/Course/fem/eigen/unsupported/Eigen/FFT"
    "F:/Course/fem/eigen/unsupported/Eigen/IterativeSolvers"
    "F:/Course/fem/eigen/unsupported/Eigen/KroneckerProduct"
    "F:/Course/fem/eigen/unsupported/Eigen/LevenbergMarquardt"
    "F:/Course/fem/eigen/unsupported/Eigen/MatrixFunctions"
    "F:/Course/fem/eigen/unsupported/Eigen/MoreVectorization"
    "F:/Course/fem/eigen/unsupported/Eigen/MPRealSupport"
    "F:/Course/fem/eigen/unsupported/Eigen/NNLS"
    "F:/Course/fem/eigen/unsupported/Eigen/NonLinearOptimization"
    "F:/Course/fem/eigen/unsupported/Eigen/NumericalDiff"
    "F:/Course/fem/eigen/unsupported/Eigen/OpenGLSupport"
    "F:/Course/fem/eigen/unsupported/Eigen/Polynomials"
    "F:/Course/fem/eigen/unsupported/Eigen/Skyline"
    "F:/Course/fem/eigen/unsupported/Eigen/SparseExtra"
    "F:/Course/fem/eigen/unsupported/Eigen/SpecialFunctions"
    "F:/Course/fem/eigen/unsupported/Eigen/Splines"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "F:/Course/fem/eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("F:/Course/fem/BUILD/debug/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

