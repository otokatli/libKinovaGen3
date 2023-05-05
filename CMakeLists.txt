cmake_minimum_required(VERSION 3.10)

project(KinovaGen3 VERSION 0.1
                   DESCRIPTION "Kinova Gen3 kinematics and dynamics"
                   LANGUAGES CXX)

set(ENV{EIGEN3_INCLUDE_DIR} /usr/include/eigen3)
find_package( Eigen3 REQUIRED )

add_library(KinovaGen3Lib STATIC src/KinovaGen3.cpp include/KinovaGen3.h)

target_include_directories(KinovaGen3Lib PUBLIC include)
target_include_directories(KinovaGen3Lib PUBLIC $ENV{EIGEN3_INCLUDE_DIR})

target_compile_definitions(KinovaGen3Lib PUBLIC cxx_std_17)

add_executable(test test/testLibrary.cpp)

target_link_libraries(test PUBLIC KinovaGen3Lib)