# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.30.0/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.30.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build"

# Include any dependencies generated for this target.
include src/komega/CMakeFiles/komega.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/komega/CMakeFiles/komega.dir/compiler_depend.make

# Include the progress variables for this target.
include src/komega/CMakeFiles/komega.dir/progress.make

# Include the compile flags for this target's objects.
include src/komega/CMakeFiles/komega.dir/flags.make

src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o: src/komega/CMakeFiles/komega.dir/flags.make
src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o: /Users/zhengbangzhou/OneDrive\ -\ University\ of\ Toronto/PhD\ Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_bicg.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_bicg.F90" -o CMakeFiles/komega.dir/komega_bicg.F90.o

src/komega/CMakeFiles/komega.dir/komega_bicg.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/komega.dir/komega_bicg.F90.i"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_bicg.F90" > CMakeFiles/komega.dir/komega_bicg.F90.i

src/komega/CMakeFiles/komega.dir/komega_bicg.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/komega.dir/komega_bicg.F90.s"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_bicg.F90" -o CMakeFiles/komega.dir/komega_bicg.F90.s

src/komega/CMakeFiles/komega.dir/komega_math.F90.o: src/komega/CMakeFiles/komega.dir/flags.make
src/komega/CMakeFiles/komega.dir/komega_math.F90.o: /Users/zhengbangzhou/OneDrive\ -\ University\ of\ Toronto/PhD\ Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_math.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object src/komega/CMakeFiles/komega.dir/komega_math.F90.o"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_math.F90" -o CMakeFiles/komega.dir/komega_math.F90.o

src/komega/CMakeFiles/komega.dir/komega_math.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/komega.dir/komega_math.F90.i"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_math.F90" > CMakeFiles/komega.dir/komega_math.F90.i

src/komega/CMakeFiles/komega.dir/komega_math.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/komega.dir/komega_math.F90.s"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_math.F90" -o CMakeFiles/komega.dir/komega_math.F90.s

src/komega/CMakeFiles/komega.dir/komega_vals.F90.o: src/komega/CMakeFiles/komega.dir/flags.make
src/komega/CMakeFiles/komega.dir/komega_vals.F90.o: /Users/zhengbangzhou/OneDrive\ -\ University\ of\ Toronto/PhD\ Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_vals.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object src/komega/CMakeFiles/komega.dir/komega_vals.F90.o"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_vals.F90" -o CMakeFiles/komega.dir/komega_vals.F90.o

src/komega/CMakeFiles/komega.dir/komega_vals.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/komega.dir/komega_vals.F90.i"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_vals.F90" > CMakeFiles/komega.dir/komega_vals.F90.i

src/komega/CMakeFiles/komega.dir/komega_vals.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/komega.dir/komega_vals.F90.s"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega/komega_vals.F90" -o CMakeFiles/komega.dir/komega_vals.F90.s

# Object files for target komega
komega_OBJECTS = \
"CMakeFiles/komega.dir/komega_bicg.F90.o" \
"CMakeFiles/komega.dir/komega_math.F90.o" \
"CMakeFiles/komega.dir/komega_vals.F90.o"

# External object files for target komega
komega_EXTERNAL_OBJECTS =

src/komega/libkomega.dylib: src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o
src/komega/libkomega.dylib: src/komega/CMakeFiles/komega.dir/komega_math.F90.o
src/komega/libkomega.dylib: src/komega/CMakeFiles/komega.dir/komega_vals.F90.o
src/komega/libkomega.dylib: src/komega/CMakeFiles/komega.dir/build.make
src/komega/libkomega.dylib: /opt/homebrew/Cellar/mpich/4.2.2/lib/libmpifort.dylib
src/komega/libkomega.dylib: /opt/homebrew/Cellar/mpich/4.2.2/lib/libmpi.dylib
src/komega/libkomega.dylib: /opt/homebrew/Cellar/mpich/4.2.2/lib/libpmpi.dylib
src/komega/libkomega.dylib: src/komega/CMakeFiles/komega.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Linking Fortran shared library libkomega.dylib"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/komega.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/komega/CMakeFiles/komega.dir/build: src/komega/libkomega.dylib
.PHONY : src/komega/CMakeFiles/komega.dir/build

src/komega/CMakeFiles/komega.dir/clean:
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" && $(CMAKE_COMMAND) -P CMakeFiles/komega.dir/cmake_clean.cmake
.PHONY : src/komega/CMakeFiles/komega.dir/clean

src/komega/CMakeFiles/komega.dir/depend:
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/src/komega" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/src/komega/CMakeFiles/komega.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : src/komega/CMakeFiles/komega.dir/depend

