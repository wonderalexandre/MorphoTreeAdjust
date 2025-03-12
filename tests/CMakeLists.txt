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


add_executable(testePruningLeaf testePruningLeaf.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustment.cpp
)

add_executable(TestComponentTree TestComponentTree.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustment.cpp
)

add_executable(dgmm2025 dgmm2025.cpp
    ../morphoTreeAdjust/src/AdjacencyRelation.cpp
    ../morphoTreeAdjust/src/ComponentTreeAdjustment.cpp
    ./external/stb/stb_image.cpp
)
target_compile_options(dgmm2025 PRIVATE -Wno-deprecated-declarations)
