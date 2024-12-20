project(Teste)

cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


include_directories(
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/morphoTreeAdjust/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/external/stb"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests"
  )


# Exemplos com a morphotree
add_executable(teste teste.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/NodeCT.cpp
    ../morphoTreeAdjust/src/ComponentTree.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustment.cpp
)

# Exemplos com a morphotree
add_executable(dgmm gdmm2024.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/NodeCT.cpp
    ../morphoTreeAdjust/src/ComponentTree.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustment.cpp
)