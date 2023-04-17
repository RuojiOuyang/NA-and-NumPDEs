CMAKE_MINIMUM_REQUIRED(VERSION 3.2) 

PROJECT(Multigrid)  

LINK_DIRECTORIES(
    "/usr/local/lib"
    "${PROJECT_SOURCE_DIR}/lib"
)
# LINK_LIBRARIES(lapack lapacke refblas)
LINK_LIBRARIES(lapack)


INCLUDE_DIRECTORIES(
    "${PROJECT_SOURCE_DIR}/include"
    "/usr/local/include"
    "/usr/local/include/eigen3"
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17") 

#default compiler
SET(CMAKE_C_COMPILER "/usr/bin/gcc") 
SET(CMAKE_CXX_COMPILER "/usr/bin/g++") 

#SET(CMAKE_BUILD_TYPE "Debug")  # Debug模式 选项: Release Debug MinSizeRel RelWithDebInfo

#SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")  # debug模式下 gdb相关选项
#SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall") # release模式下 gdb相关选项

# set(CMAKE_VERBOSE_MAKEFILE ON) # 开启调试 出现问题时开启

#enable_testing() # 打开测试

add_definitions(-DDEBUG) # 相当于代码中 #define DEBUG

#add_subdirectory(test) # 添加test子目录

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin") # 可执行文件输出目录
SET(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib") 

add_custom_target(doc
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND doxygen Doxyfile
)

set (Header "${PROJECT_SOURCE_DIR}/include")
add_library(genJson
    SHARED
    ${Header}/genJson.h
    ${Header}/genJson.cpp
    ${Header}/global.cpp
    ${Header}/global.h
)
target_link_libraries(genJson
    -ljsoncpp
)

add_library(Multigrid
    SHARED
    ${Header}/global.cpp
    ${Header}/global.h
    ${Header}/Vec2d.h
    ${Header}/conditions.h
    ${Header}/conditions.cpp
    ${Header}/Generator.h
    ${Header}/Generator.cpp
    ${Header}/operator.h
    ${Header}/operator.cpp
    ${Header}/Multigrid.h
    ${Header}/Multigrid.cpp
    # ${Header}/MGsolver.h
    # ${Header}/MGsolver.cpp
)
target_link_libraries(Multigrid
    -ljsoncpp
)

add_executable(test1
    ./include/test_helper.h
    ./include/test_helper.cpp
    ./src/test1.cpp
)
target_link_libraries(test1
    -ljsoncpp
    -lgenJson
    -lMultigrid
)

add_executable(test2
    ./include/test_helper.h
    ./include/test_helper.cpp
    ./src/test2.cpp
)
target_link_libraries(test2
    -ljsoncpp
    -lgenJson
    -lMultigrid
)

