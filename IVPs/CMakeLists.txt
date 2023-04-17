CMAKE_MINIMUM_REQUIRED(VERSION 3.2) 

PROJECT(IVPs)  

LINK_DIRECTORIES(
    "/usr/local/lib"
    "${PROJECT_SOURCE_DIR}/lib"
)
# LINK_LIBRARIES(lapack lapacke refblas)
LINK_LIBRARIES(lapack)

INCLUDE_DIRECTORIES(
    "${PROJECT_SOURCE_DIR}/include"
    "/usr/local/include"
    "/usr/local/include/eigen3/unsupported"
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



aux_source_directory ("${PROJECT_SOURCE_DIR}/include" Header)
aux_source_directory ("${PROJECT_SOURCE_DIR}/include/RHSFunction" Header2)
aux_source_directory ("${PROJECT_SOURCE_DIR}/include/TimeIntegrator" Header3)
aux_source_directory ("${PROJECT_SOURCE_DIR}/include/TimeIntegrator/RungeKutta" Header4)
aux_source_directory ("${PROJECT_SOURCE_DIR}/include/TimeIntegrator/LinearMultistep" Header5)

# add_executable(mytest 
#    ./src/test.cpp
# )
# target_link_libraries(mytest
#         -ljsoncpp
#         ${PROJECT_SOURCE_DIR}/build/libFD.so
# )



add_library(IVPs
    SHARED
    ${Header}
    ${Header2}
    ${Header3}
    ${Header4}
    ${Header5}
)
target_link_libraries(IVPs
    -ljsoncpp
)

add_executable(cr_test
    ./src/cr_test.cpp
)
target_link_libraries(cr_test
    -lIVPs
    -ljsoncpp
)

add_executable(main
    ./src/main.cpp
)
target_link_libraries(main
    -lIVPs
    -ljsoncpp
)