# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Note that incremental build could trigger a call to cmake_copy_f90_mod on each re-build
src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o: src/komega/CMakeFiles/komega.dir/komega_math.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o: src/komega/CMakeFiles/komega.dir/komega_parameter.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o: src/komega/CMakeFiles/komega.dir/komega_vals_c.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o: src/komega/CMakeFiles/komega.dir/komega_vecs_c.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o.provides.build: src/komega/CMakeFiles/komega.dir/komega_bicg.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_bicg.mod.stamp: src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod src/komega/komega_bicg.mod src/komega/CMakeFiles/komega.dir/komega_bicg.mod.stamp GNU
src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o.provides.build
src/komega/CMakeFiles/komega.dir/build: src/komega/CMakeFiles/komega.dir/komega_bicg.F90.o.provides.build
src/komega/CMakeFiles/komega.dir/komega_math.F90.o: src/komega/CMakeFiles/komega.dir/komega_parameter.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_math.F90.o: /opt/homebrew/Cellar/mpich/4.2.2/include/mpi.mod
src/komega/CMakeFiles/komega.dir/komega_math.F90.o.provides.build: src/komega/CMakeFiles/komega.dir/komega_math.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_math.mod.stamp: src/komega/CMakeFiles/komega.dir/komega_math.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod src/komega/komega_math.mod src/komega/CMakeFiles/komega.dir/komega_math.mod.stamp GNU
src/komega/CMakeFiles/komega.dir/komega_math.F90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/komega/CMakeFiles/komega.dir/komega_math.F90.o.provides.build
src/komega/CMakeFiles/komega.dir/build: src/komega/CMakeFiles/komega.dir/komega_math.F90.o.provides.build
src/komega/CMakeFiles/komega.dir/komega_vals.F90.o.provides.build: src/komega/CMakeFiles/komega.dir/komega_parameter.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_parameter.mod.stamp: src/komega/CMakeFiles/komega.dir/komega_vals.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod src/komega/komega_parameter.mod src/komega/CMakeFiles/komega.dir/komega_parameter.mod.stamp GNU
src/komega/CMakeFiles/komega.dir/komega_vals.F90.o.provides.build: src/komega/CMakeFiles/komega.dir/komega_vals_c.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_vals_c.mod.stamp: src/komega/CMakeFiles/komega.dir/komega_vals.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod src/komega/komega_vals_c.mod src/komega/CMakeFiles/komega.dir/komega_vals_c.mod.stamp GNU
src/komega/CMakeFiles/komega.dir/komega_vals.F90.o.provides.build: src/komega/CMakeFiles/komega.dir/komega_vals_r.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_vals_r.mod.stamp: src/komega/CMakeFiles/komega.dir/komega_vals.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod src/komega/komega_vals_r.mod src/komega/CMakeFiles/komega.dir/komega_vals_r.mod.stamp GNU
src/komega/CMakeFiles/komega.dir/komega_vals.F90.o.provides.build: src/komega/CMakeFiles/komega.dir/komega_vecs_c.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_vecs_c.mod.stamp: src/komega/CMakeFiles/komega.dir/komega_vals.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod src/komega/komega_vecs_c.mod src/komega/CMakeFiles/komega.dir/komega_vecs_c.mod.stamp GNU
src/komega/CMakeFiles/komega.dir/komega_vals.F90.o.provides.build: src/komega/CMakeFiles/komega.dir/komega_vecs_r.mod.stamp
src/komega/CMakeFiles/komega.dir/komega_vecs_r.mod.stamp: src/komega/CMakeFiles/komega.dir/komega_vals.F90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod src/komega/komega_vecs_r.mod src/komega/CMakeFiles/komega.dir/komega_vecs_r.mod.stamp GNU
src/komega/CMakeFiles/komega.dir/komega_vals.F90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/komega/CMakeFiles/komega.dir/komega_vals.F90.o.provides.build
src/komega/CMakeFiles/komega.dir/build: src/komega/CMakeFiles/komega.dir/komega_vals.F90.o.provides.build
