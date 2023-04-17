CMAKE_MINIMUM_REQUIRED(VERSION 3.2) 

PROJECT(FD_Method)  

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

# add_custom_target(plot
#     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/octave
#     COMMAND octave problem_B.m 
#     COMMAND octave problem_C.m 
#     COMMAND octave problem_E.m 
#     COMMAND octave problem_F.m 
# )

# add_custom_target(clean_all
#     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#     COMMAND rm -R ./octave/*.m ./octave/images/* ./build/* ./bin/*
# )

# add_custom_target(clean_images
#     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/octave/images
#     COMMAND rm *
# )

# add_custom_target(exec_all
#     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
#     COMMAND ./problem_B
#     COMMAND ./problem_CandD
#     COMMAND ./problem_E
#     COMMAND ./problem_F
# )

add_custom_target(doc
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND doxygen Doxyfile
)

add_custom_target(genJson
    COMMAND ../bin/myJson1
    COMMAND ../bin/myJson2
    COMMAND ../bin/myJson3
    COMMAND ../bin/myJson4
    COMMAND ../bin/myJson5
    COMMAND ../bin/myJson6
)

add_custom_target(start_test
    COMMAND ../bin/myJson1
    COMMAND ../bin/myJson2
    COMMAND ../bin/myJson3
    COMMAND ../bin/myJson4
    COMMAND ../bin/myJson5
    COMMAND ../bin/myJson6
    COMMAND ../bin/myJson7
    COMMAND ../bin/myJson8
    COMMAND ../bin/myJson9
    COMMAND ../bin/mytest
)

FILE(GLOB_RECURSE mains RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/myJson*.cpp")
# aux_source_directory ("${PROJECT_SOURCE_DIR}/include" Header)
set (Header "${PROJECT_SOURCE_DIR}/include")
add_library(FD
    SHARED
    ${Header}/conditions.cpp  
    ${Header}/conditions.h
    ${Header}/FD.cpp  
    ${Header}/FD.h
    ${Header}/global.cpp
    ${Header}/global.h
    ${Header}/Vec2d.h
    ${Header}/genJson.h
    ${Header}/genJson.cpp
)
target_link_libraries(FD
    -ljsoncpp
)


FOREACH(mainfile IN LISTS mains)
    # Get file name without directory
    get_filename_component(mainname ${mainfile} NAME_WE)
    add_executable(${mainname} 
        ${mainfile}
        # ${Header}
    )
    target_link_libraries(${mainname}  
        -ljsoncpp
        ${PROJECT_SOURCE_DIR}/build/libFD.so
    )
ENDFOREACH()

# aux_source_directory ("${PROJECT_SOURCE_DIR}/include" Header)


add_executable(mytest 
   ./src/test.cpp
)
target_link_libraries(mytest
        -ljsoncpp
        ${PROJECT_SOURCE_DIR}/build/libFD.so
)

# add_executable(main 
#    ./src/main.cpp
# )
# target_link_libraries(main
#         -ljsoncpp
#         ${PROJECT_SOURCE_DIR}/build/libFD.so
# )

# add_executable(main 
#    ./src/test.cpp
#    ${Header}
# )
# target_link_libraries(main
#         -ljsoncpp
#         -llapacke 
#         -llapack 
#         -lrefblas
#         -lgfortran
# )

