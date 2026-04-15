# build
# cmake -DPYTHON_LIBRARY_DIR=$(python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))") -DPYTHON_EXECUTABLE=$(python -c "import sys; print(sys.executable)") ../
# make
# make install
cmake_minimum_required(VERSION 3.10)
include("${CMAKE_CURRENT_LIST_DIR}/cmake/AppleSdkFix.cmake")
project(MorphoTreeAdjust LANGUAGES CXX)
include(GNUInstallDirs)

# C++20 without compiler-specific extensions
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Define the default build type if none is specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# Build-type-specific flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wextra -Wpedantic")

# Helpful status messages
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")

option(MTA_INSTALL_PUBLIC_CORE_HEADERS "Install only the public C++ core headers." ON)
set(MTA_PYTHON_PACKAGE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/morphoTreeAdjust" CACHE PATH
  "Install directory for the Python package relative to the install prefix.")

include_directories("${CMAKE_SOURCE_DIR}/python")

file (GLOB SOURCE_FILES "morphoTreeAdjust/src/*.cpp")
file (GLOB HEADER_FILES "morphoTreeAdjust/include/*.hpp" "morphoTreeAdjust/include/*.tpp")
set(PYTHON_FILES "morphoTreeAdjust/morphoTreeAdjust.cpp")
set(MTA_PYTHON_PACKAGE_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/python/morphoTreeAdjust/__init__.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/python/morphoTreeAdjust/utils.py"
)

set(MTA_PUBLIC_CORE_HEADER_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include/MorphoTreeAdjust.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include/AdjacencyRelation.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include/AttributeComputer.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include/Common.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include/ComponentTreeCasf.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include/DynamicComponentTree.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include/DynamicComponentTreeAdjustment.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include/DynamicComponentTreeAdjustmentLeaf.hpp"
)

# Set up such that XCode organizes the files
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )

include(pybind11.cmake)
pybind11_add_module(morphoTreeAdjust 
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)

target_include_directories(morphoTreeAdjust PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/MorphoTreeAdjust>
)

target_link_libraries(morphoTreeAdjust PUBLIC)

message(STATUS
  "MorphoTreeAdjust install surface:"
  " python-package-dir=${MTA_PYTHON_PACKAGE_INSTALL_DIR}"
  " public-core-headers=${MTA_INSTALL_PUBLIC_CORE_HEADERS}"
  " transitional-headers=OFF"
)

install(TARGETS morphoTreeAdjust
  COMPONENT python
  LIBRARY DESTINATION "${MTA_PYTHON_PACKAGE_INSTALL_DIR}"
  )

install(FILES ${MTA_PYTHON_PACKAGE_FILES}
  COMPONENT python
  DESTINATION "${MTA_PYTHON_PACKAGE_INSTALL_DIR}"
)

if(MTA_INSTALL_PUBLIC_CORE_HEADERS)
  install(FILES ${MTA_PUBLIC_CORE_HEADER_FILES}
    COMPONENT core_headers
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/MorphoTreeAdjust"
  )
endif()
