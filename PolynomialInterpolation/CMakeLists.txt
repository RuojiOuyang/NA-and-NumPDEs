CMAKE_MINIMUM_REQUIRED(VERSION 3.2) 

PROJECT(Polynomail_Interpolation)  

LINK_DIRECTORIES("/usr/local/lib")
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") # 添加c++11标准支持

#default compiler
SET(CMAKE_C_COMPILER "/usr/bin/gcc") 
SET(CMAKE_CXX_COMPILER "/usr/bin/g++") 

#SET(CMAKE_BUILD_TYPE "Debug")  # Debug模式 选项: Release Debug MinSizeRel RelWithDebInfo

#SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")  # debug模式下 gdb相关选项
#SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall") # release模式下 gdb相关选项

# set(CMAKE_VERBOSE_MAKEFILE ON) # 开启调试 出现问题时开启

enable_testing() # 打开测试

add_definitions(-DDEBUG) # 相当于代码中 #define DEBUG

add_subdirectory(test) # 添加test子目录

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin") # 可执行文件输出目录

add_custom_target(plot
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/plot
    COMMAND octave problem_B.m 
    COMMAND octave problem_C.m 
    COMMAND octave problem_D.m 
    COMMAND octave problem_E.m 
)

add_custom_target(clean_output
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/output
    COMMAND rm ./images/*.png
    COMMAND rm ./4octave/*.txt
)

FILE(GLOB_RECURSE mains RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/problem_*.cpp")
FOREACH(mainfile IN LISTS mains)
    # Get file name without directory
    get_filename_component(mainname ${mainfile} NAME_WE)
    add_executable(${mainname} 
        ${mainfile}
        ./src/PolyInterp.cpp
    )
ENDFOREACH()
# 也可以这么写
# add_executable(chess  # 输出名为chess的可执行文件
#    ./src/main.cpp
#    ./src/utils.cpp
# )
