CMAKE_MINIMUM_REQUIRED(VERSION 3.2) 

PROJECT(Approximation)  

LINK_DIRECTORIES("/usr/local/lib")
LINK_LIBRARIES(lapack lapacke refblas)

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include")

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

add_custom_target(doc
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND doxygen Doxyfile
    #COMMAND ln -s ./html/index.html ./doc/doc.html
)

FILE(GLOB_RECURSE mains RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/problem_*.cpp")
FOREACH(mainfile IN LISTS mains)
    # Get file name without directory
    get_filename_component(mainname ${mainfile} NAME_WE)
    add_executable(${mainname} 
        ${mainfile}
        ./src/header/DLS.cpp
    )
    target_link_libraries(${mainname}  
        -llapacke 
        -llapack 
        -lrefblas
        -lgfortran
    )
ENDFOREACH()


# add_executable(main 
#    ./src/test.cpp
# )
# target_link_libraries(main
#         -llapacke 
#         -llapack 
#         -lrefblas
#         -lgfortran
# )

