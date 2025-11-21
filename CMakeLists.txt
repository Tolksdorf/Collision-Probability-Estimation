# Written by Leon Tolksdorf, 2025. Please refer to: 

# @article{CollisionProbabilityEstimation,
#   title   = {Collision Probability Estimation for Optimization-based Vehicular Motion Planning},
#   author  = {Tolksdorf, Leon and Tejada, Arturo and Birkner, Christian and van de Wouw, Nathan},
#   journal = {arXiv preprint arXiv:2505.21161} ,
#   year    = {2025}
# }

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(_poc_function)

find_package(casadi REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)

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
    casadi
    Eigen3::Eigen)

target_include_directories(${PROJECT_NAME} PRIVATE include)

add_executable(testbench.out testbench.cpp)

target_link_libraries(testbench.out ${PROJECT_NAME} casadi Eigen3::Eigen)

set_target_properties(testbench.out PROPERTIES 
    BUILD_RPATH $ORIGIN
    INSTALL_RPATH $ORIGIN
)