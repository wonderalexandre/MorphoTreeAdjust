
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


# Diretórios de includes
include_directories(
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/external/stb"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests"
  )

# Exemplos com a morphotree
add_executable(testeRemovingAnyNode testeRemovingAnyNode.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustmentByFlatzone.cpp
)

add_executable(dgmm2025 dgmm2025.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustmentByLeaf.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustmentBySubtree.cpp
    ./external/stb/stb_image.cpp
)
target_compile_options(dgmm2025 PRIVATE -Wno-deprecated-declarations)


#
# Exemplos com a morphotree
add_executable(testePruningLeaf testePruningLeaf.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustmentByLeaf.cpp
    ./external/stb/stb_image.cpp
)

add_executable(TestComponentTree TestComponentTree.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustmentByLeaf.cpp
)

# Exemplos com a morphotree
add_executable(testePruningSubtree testePruningSubtree.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustmentBySubtree.cpp
)

#add_executable(testeComponentTreeGeneralAdjustment testeComponentTreeGeneralAdjustment.cpp
#    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
#    ../morphoTreeAdjust/src/ComponentTreeAdjustment.cpp
#)

