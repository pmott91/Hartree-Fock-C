cmake_minimum_required(VERSION 3.15.2)
project (parquant)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_STANDARD_REQUIRED_ON)

#Bring in the header files.
include_directories(include
			${EIGEN_INCLUDE_PATH}
                     ${LIBINT_INCLUDE_PATH}
			)

#Add the source files
file(GLOB SOURCES "src/*.cpp")

find_path(EIGEN_INCLUDE_PATH | NAMES Eigen | HINTS ${EIGEN_ROOT} NO_DEFAULT_PATH)
message("Eigen path is ${EIGEN_INCLUDE_PATH}")
if(EXISTS ${EIGEN_INCLUDE_PATH})
	message("-- Detecting Eigen ${EIGEN_INCLUDE_PATH} - done")
else()
	message(FATAL_ERROR "Eigen not found! Please select a proper location for the Eigen installation.")
endif(EXISTS ${EIGEN_INCLUDE_PATH})

find_path(LIBINT_INCLUDE_PATH | NAMES libint2 | HINTS ${LIBINT_ROOT} NO_DEFAULT_PATH)
message("Eigen path is ${LIBINT_INCLUDE_PATH}")
if(EXISTS ${LIBINT_INCLUDE_PATH})
       message("-- Detecting libint ${LIBINT_INCLUDE_PATH} - done")
else()
       message(FATAL_ERROR "Libint not found! Please select a proper location for the Libint installation.")
endif(EXISTS ${LIBINT_INCLUDE_PATH})
#find_package (Eigen3 3.3 NO_MODULE)

add_executable(HF ${SOURCES}) 

#target_link_libraries(HF ${EIGEN_INCLUDE_PATH}::Eigen)
