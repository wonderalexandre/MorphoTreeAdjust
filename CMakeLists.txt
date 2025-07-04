#compilar
# cmake -DPYTHON_LIBRARY_DIR=$(python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))") -DPYTHON_EXECUTABLE=$(python -c "import sys; print(sys.executable)") ../
# make
# make install
cmake_minimum_required(VERSION 3.10)
project(MorphoTreeAdjust LANGUAGES CXX)

# C++20 sem extensões específicas de compilador
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Define build padrão, se não for especificado
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# Flags específicas por tipo de build
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra -Wpedantic")

# Mensagens úteis
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")



include_directories("${CMAKE_SOURCE_DIR}/python")

file (GLOB SOURCE_FILES "morphoTreeAdjust/src/*.cpp" "morphoTreeAdjust/pybind/*.cpp")
file (GLOB HEADER_FILES "morphoTreeAdjust/include/*.hpp" "morphoTreeAdjust/include/*.tpp" "morphoTreeAdjust/pybind/*.hpp")
file (GLOB PYTHON_FILES "morphoTreeAdjust/*.cpp" "morphoTreeAdjust/*.hpp")

# Set up such that XCode organizes the files
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )

include(pybind11.cmake)
pybind11_add_module(morphoTreeAdjust 
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)

target_link_libraries(morphoTreeAdjust PUBLIC)

install(TARGETS morphoTreeAdjust
  COMPONENT python
  LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
  )
