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
include tool/CMakeFiles/greenr2k.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tool/CMakeFiles/greenr2k.dir/compiler_depend.make

# Include the progress variables for this target.
include tool/CMakeFiles/greenr2k.dir/progress.make

# Include the compile flags for this target's objects.
include tool/CMakeFiles/greenr2k.dir/flags.make

tool/CMakeFiles/greenr2k.dir/greenr2k.F90.o: tool/CMakeFiles/greenr2k.dir/flags.make
tool/CMakeFiles/greenr2k.dir/greenr2k.F90.o: /Users/zhengbangzhou/OneDrive\ -\ University\ of\ Toronto/PhD\ Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/tool/greenr2k.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object tool/CMakeFiles/greenr2k.dir/greenr2k.F90.o"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/tool" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/tool/greenr2k.F90" -o CMakeFiles/greenr2k.dir/greenr2k.F90.o

tool/CMakeFiles/greenr2k.dir/greenr2k.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/greenr2k.dir/greenr2k.F90.i"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/tool" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/tool/greenr2k.F90" > CMakeFiles/greenr2k.dir/greenr2k.F90.i

tool/CMakeFiles/greenr2k.dir/greenr2k.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/greenr2k.dir/greenr2k.F90.s"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/tool" && /opt/homebrew/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/tool/greenr2k.F90" -o CMakeFiles/greenr2k.dir/greenr2k.F90.s

# Object files for target greenr2k
greenr2k_OBJECTS = \
"CMakeFiles/greenr2k.dir/greenr2k.F90.o"

# External object files for target greenr2k
greenr2k_EXTERNAL_OBJECTS =

tool/greenr2k: tool/CMakeFiles/greenr2k.dir/greenr2k.F90.o
tool/greenr2k: tool/CMakeFiles/greenr2k.dir/build.make
tool/greenr2k: tool/CMakeFiles/greenr2k.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable greenr2k"
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/tool" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/greenr2k.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tool/CMakeFiles/greenr2k.dir/build: tool/greenr2k
.PHONY : tool/CMakeFiles/greenr2k.dir/build

tool/CMakeFiles/greenr2k.dir/clean:
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/tool" && $(CMAKE_COMMAND) -P CMakeFiles/greenr2k.dir/cmake_clean.cmake
.PHONY : tool/CMakeFiles/greenr2k.dir/clean

tool/CMakeFiles/greenr2k.dir/depend:
	cd "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi-3.5.2/tool" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/tool" "/Users/zhengbangzhou/OneDrive - University of Toronto/PhD Stuff/Projects/field_induced_phase_transition/HPhi.build/tool/CMakeFiles/greenr2k.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : tool/CMakeFiles/greenr2k.dir/depend

