project(Teste)

cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


include_directories(
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/teste"
  )


# Exemplos com a morphotree
add_executable(teste teste.cpp
    ../src/AdjacencyRelation.cpp
    ../src/NodeCT.cpp
    ../src/ComponentTree.cpp
    ../src/ComponentTreeAdjustment.cpp
)
