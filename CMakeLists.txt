CMAKE_MINIMUM_REQUIRED(VERSION 2.8)  
PROJECT(circle)  
set(CMAKE_CXX_STANDARD 11)
INCLUDE_DIRECTORIES("/usr/include/eigen3") #("${EIGEN3_INCLUDE_DIR}")  
#SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")  
#IF(NOT EIGEN3_INCLUDE_DIR)  
#    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")  
#ENDIF()  
#INCLUDE_DIRECTORIES(/usr/include/eigen3) #("${EIGEN3_INCLUDE_DIR}")  
  
#ADD_EXECUTABLE(test test.cpp)  
#include_directories(src)
# 找到src目录及其子目录下的所有c文件和h文件
#file(GLOB_RECURSE PROJECT_SOURCES "src/*.cpp" "src/*.c" "src/*.cc")
set(sources src/eigen_polyfit.cpp)

add_executable(polyfit ${sources})
#file(GLOB_RECURSE PROJECT_HEADERS "src/*.hpp" "src/*.h" "src/*.hh")

#add_executable(polyfit ${PROJECT_SOURCES}  ${PROJECT_HEADERS})

target_link_libraries(polyfit python2.7)
