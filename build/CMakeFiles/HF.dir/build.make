# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /autofs/nccs-svm1_sw/summit/.swci/0-core/opt/spack/20180914/linux-rhel7-ppc64le/gcc-4.8.5/cmake-3.15.2-xit2o3iepxvqbyku77lwcugufilztu7t/bin/cmake

# The command to remove a file.
RM = /autofs/nccs-svm1_sw/summit/.swci/0-core/opt/spack/20180914/linux-rhel7-ppc64le/gcc-4.8.5/cmake-3.15.2-xit2o3iepxvqbyku77lwcugufilztu7t/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /ccs/home/plot/parquant/Parallel-Quantum-Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /ccs/home/plot/parquant/Parallel-Quantum-Code/build

# Include any dependencies generated for this target.
include CMakeFiles/HF.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/HF.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/HF.dir/flags.make

CMakeFiles/HF.dir/src/HF.cpp.o: CMakeFiles/HF.dir/flags.make
CMakeFiles/HF.dir/src/HF.cpp.o: ../src/HF.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/ccs/home/plot/parquant/Parallel-Quantum-Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/HF.dir/src/HF.cpp.o"
	/sw/summit/gcc/6.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/HF.dir/src/HF.cpp.o -c /ccs/home/plot/parquant/Parallel-Quantum-Code/src/HF.cpp

CMakeFiles/HF.dir/src/HF.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HF.dir/src/HF.cpp.i"
	/sw/summit/gcc/6.4.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /ccs/home/plot/parquant/Parallel-Quantum-Code/src/HF.cpp > CMakeFiles/HF.dir/src/HF.cpp.i

CMakeFiles/HF.dir/src/HF.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HF.dir/src/HF.cpp.s"
	/sw/summit/gcc/6.4.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /ccs/home/plot/parquant/Parallel-Quantum-Code/src/HF.cpp -o CMakeFiles/HF.dir/src/HF.cpp.s

CMakeFiles/HF.dir/src/func.cpp.o: CMakeFiles/HF.dir/flags.make
CMakeFiles/HF.dir/src/func.cpp.o: ../src/func.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/ccs/home/plot/parquant/Parallel-Quantum-Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/HF.dir/src/func.cpp.o"
	/sw/summit/gcc/6.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/HF.dir/src/func.cpp.o -c /ccs/home/plot/parquant/Parallel-Quantum-Code/src/func.cpp

CMakeFiles/HF.dir/src/func.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HF.dir/src/func.cpp.i"
	/sw/summit/gcc/6.4.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /ccs/home/plot/parquant/Parallel-Quantum-Code/src/func.cpp > CMakeFiles/HF.dir/src/func.cpp.i

CMakeFiles/HF.dir/src/func.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HF.dir/src/func.cpp.s"
	/sw/summit/gcc/6.4.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /ccs/home/plot/parquant/Parallel-Quantum-Code/src/func.cpp -o CMakeFiles/HF.dir/src/func.cpp.s

CMakeFiles/HF.dir/src/quantnum.cpp.o: CMakeFiles/HF.dir/flags.make
CMakeFiles/HF.dir/src/quantnum.cpp.o: ../src/quantnum.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/ccs/home/plot/parquant/Parallel-Quantum-Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/HF.dir/src/quantnum.cpp.o"
	/sw/summit/gcc/6.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/HF.dir/src/quantnum.cpp.o -c /ccs/home/plot/parquant/Parallel-Quantum-Code/src/quantnum.cpp

CMakeFiles/HF.dir/src/quantnum.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HF.dir/src/quantnum.cpp.i"
	/sw/summit/gcc/6.4.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /ccs/home/plot/parquant/Parallel-Quantum-Code/src/quantnum.cpp > CMakeFiles/HF.dir/src/quantnum.cpp.i

CMakeFiles/HF.dir/src/quantnum.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HF.dir/src/quantnum.cpp.s"
	/sw/summit/gcc/6.4.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /ccs/home/plot/parquant/Parallel-Quantum-Code/src/quantnum.cpp -o CMakeFiles/HF.dir/src/quantnum.cpp.s

# Object files for target HF
HF_OBJECTS = \
"CMakeFiles/HF.dir/src/HF.cpp.o" \
"CMakeFiles/HF.dir/src/func.cpp.o" \
"CMakeFiles/HF.dir/src/quantnum.cpp.o"

# External object files for target HF
HF_EXTERNAL_OBJECTS =

HF: CMakeFiles/HF.dir/src/HF.cpp.o
HF: CMakeFiles/HF.dir/src/func.cpp.o
HF: CMakeFiles/HF.dir/src/quantnum.cpp.o
HF: CMakeFiles/HF.dir/build.make
HF: CMakeFiles/HF.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/ccs/home/plot/parquant/Parallel-Quantum-Code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable HF"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/HF.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/HF.dir/build: HF

.PHONY : CMakeFiles/HF.dir/build

CMakeFiles/HF.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/HF.dir/cmake_clean.cmake
.PHONY : CMakeFiles/HF.dir/clean

CMakeFiles/HF.dir/depend:
	cd /ccs/home/plot/parquant/Parallel-Quantum-Code/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ccs/home/plot/parquant/Parallel-Quantum-Code /ccs/home/plot/parquant/Parallel-Quantum-Code /ccs/home/plot/parquant/Parallel-Quantum-Code/build /ccs/home/plot/parquant/Parallel-Quantum-Code/build /ccs/home/plot/parquant/Parallel-Quantum-Code/build/CMakeFiles/HF.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/HF.dir/depend

