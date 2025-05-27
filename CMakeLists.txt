cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
# set(CMAKE_C_COMPILER "gcc-11.4.0")
# set(CMAKE_CXX_COMPILER "/usr/bin/g++-7.5")

project(_poc_function)

set(${PROJECT_NAME}.out CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(${PROJECT_NAME} SHARED 
    collision_angle_intervals.cpp
    side_functions.cpp
    int_unions.cpp
    poc_overhead_fu.cpp
    sort_intervals.cpp
    twoD_integral.cpp
    MCS_rectangle.cpp
)
include_directories(include)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
    casadi)

target_include_directories(${PROJECT_NAME} PRIVATE include)

add_executable(testbench.out testbench.cpp)

target_link_libraries(testbench.out ${PROJECT_NAME})

#uncomment if the testbenach does have additional header files
#target_include_directories(testbench.out PRIVATE include)

set_target_properties(testbench.out PROPERTIES 
    BUILD_RPATH $ORIGIN
    INSTALL_RPATH $ORIGIN
)