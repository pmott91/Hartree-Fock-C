cmake_minimum_required(VERSION 3.16.1)
project (parquant)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_STANDARD_REQUIRED_ON)

#Bring in the header files.
include_directories(include)

#Add the source files
file(GLOB SOURCES "src/*.cpp")

find_package (Eigen3 3.3 NO_MODULE)

add_executable(HF ${SOURCES}) 

target_link_libraries(HF Eigen3::Eigen)
