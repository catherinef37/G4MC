# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/Nickie/JLab/HallA/G4MC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/Nickie/JLab/HallA/G4MC/cmake

# Include any dependencies generated for this target.
include CMakeFiles/XSModel.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/XSModel.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/XSModel.dir/flags.make

CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o: ../XSModel/Bosted/PBosted.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o -c /home/Nickie/JLab/HallA/G4MC/XSModel/Bosted/PBosted.cc

CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Nickie/JLab/HallA/G4MC/XSModel/Bosted/PBosted.cc > CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.i

CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Nickie/JLab/HallA/G4MC/XSModel/Bosted/PBosted.cc -o CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.s

CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o.requires

CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o.provides: CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o.provides

CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o

CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o: ../XSModel/Bosted/bosted.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/Nickie/JLab/HallA/G4MC/XSModel/Bosted/bosted.f -o CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o

CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o.requires

CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o.provides: CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o.provides

CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o

CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o: ../XSModel/Compton/RCSXS.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o -c /home/Nickie/JLab/HallA/G4MC/XSModel/Compton/RCSXS.cc

CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Nickie/JLab/HallA/G4MC/XSModel/Compton/RCSXS.cc > CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.i

CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Nickie/JLab/HallA/G4MC/XSModel/Compton/RCSXS.cc -o CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.s

CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o.requires

CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o.provides: CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o.provides

CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o

CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o: ../XSModel/Elas/ElasModel.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o -c /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/ElasModel.cc

CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/ElasModel.cc > CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.i

CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/ElasModel.cc -o CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.s

CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o.requires

CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o.provides: CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o.provides

CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o

CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o: ../XSModel/Elas/gauss_legendre.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o -c /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/gauss_legendre.cc

CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/gauss_legendre.cc > CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.i

CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/gauss_legendre.cc -o CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.s

CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o.requires

CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o.provides: CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o.provides

CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o

CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o: ../XSModel/Elas/deq2.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/deq2.f -o CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o

CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o.requires

CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o.provides: CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o.provides

CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o

CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o: ../XSModel/Elas/carbon_larry.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/carbon_larry.f -o CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o

CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o.requires

CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o.provides: CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o.provides

CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o

CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o: ../XSModel/Elas/chgden-fb.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/chgden-fb.f -o CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o

CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o.requires

CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o.provides: CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o.provides

CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o

CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o: ../XSModel/Elas/indeq2.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/indeq2.f -o CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o

CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o.requires

CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o.provides: CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o.provides

CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o

CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o: ../XSModel/Elas/ravpro.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/Nickie/JLab/HallA/G4MC/XSModel/Elas/ravpro.f -o CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o

CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o.requires

CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o.provides: CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o.provides

CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o: ../XSModel/QFS_N_EPC/QFS_N_EPC.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_11)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o -c /home/Nickie/JLab/HallA/G4MC/XSModel/QFS_N_EPC/QFS_N_EPC.cc

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Nickie/JLab/HallA/G4MC/XSModel/QFS_N_EPC/QFS_N_EPC.cc > CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.i

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Nickie/JLab/HallA/G4MC/XSModel/QFS_N_EPC/QFS_N_EPC.cc -o CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.s

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o.requires

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o.provides: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o.provides

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o.provides.build: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o: ../XSModel/QFS_N_EPC/EPC_photon.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_12)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o -c /home/Nickie/JLab/HallA/G4MC/XSModel/QFS_N_EPC/EPC_photon.cc

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Nickie/JLab/HallA/G4MC/XSModel/QFS_N_EPC/EPC_photon.cc > CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.i

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Nickie/JLab/HallA/G4MC/XSModel/QFS_N_EPC/EPC_photon.cc -o CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.s

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o.requires

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o.provides: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o.provides

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o.provides.build: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o: ../XSModel/QFS_N_EPC/epc.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_13)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/Nickie/JLab/HallA/G4MC/XSModel/QFS_N_EPC/epc.f -o CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o.requires

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o.provides: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o.provides

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o.provides.build: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o: ../XSModel/QFS_N_EPC/qfs.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_14)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/Nickie/JLab/HallA/G4MC/XSModel/QFS_N_EPC/qfs.f -o CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o.requires

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o.provides: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o.provides

CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o.provides.build: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o

CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o: ../XSModel/Wiser/WISER_c.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_15)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o -c /home/Nickie/JLab/HallA/G4MC/XSModel/Wiser/WISER_c.cc

CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/Nickie/JLab/HallA/G4MC/XSModel/Wiser/WISER_c.cc > CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.i

CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/Nickie/JLab/HallA/G4MC/XSModel/Wiser/WISER_c.cc -o CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.s

CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o.requires

CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o.provides: CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o.provides

CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o

CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o: CMakeFiles/XSModel.dir/flags.make
CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o: ../XSModel/Wiser/wiser.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles $(CMAKE_PROGRESS_16)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/Nickie/JLab/HallA/G4MC/XSModel/Wiser/wiser.f -o CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o

CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o.requires:
.PHONY : CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o.requires

CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o.provides: CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o.requires
	$(MAKE) -f CMakeFiles/XSModel.dir/build.make CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o.provides.build
.PHONY : CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o.provides

CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o.provides.build: CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o

# Object files for target XSModel
XSModel_OBJECTS = \
"CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o" \
"CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o" \
"CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o" \
"CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o" \
"CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o" \
"CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o" \
"CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o" \
"CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o" \
"CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o" \
"CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o" \
"CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o" \
"CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o" \
"CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o" \
"CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o" \
"CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o" \
"CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o"

# External object files for target XSModel
XSModel_EXTERNAL_OBJECTS =

libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o
libXSModel.a: CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o
libXSModel.a: CMakeFiles/XSModel.dir/build.make
libXSModel.a: CMakeFiles/XSModel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libXSModel.a"
	$(CMAKE_COMMAND) -P CMakeFiles/XSModel.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/XSModel.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/XSModel.dir/build: libXSModel.a
.PHONY : CMakeFiles/XSModel.dir/build

CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Bosted/PBosted.cc.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Bosted/bosted.f.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Compton/RCSXS.cc.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Elas/ElasModel.cc.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Elas/gauss_legendre.cc.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Elas/deq2.f.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Elas/carbon_larry.f.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Elas/chgden-fb.f.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Elas/indeq2.f.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Elas/ravpro.f.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/QFS_N_EPC.cc.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/EPC_photon.cc.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/epc.f.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/QFS_N_EPC/qfs.f.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Wiser/WISER_c.cc.o.requires
CMakeFiles/XSModel.dir/requires: CMakeFiles/XSModel.dir/XSModel/Wiser/wiser.f.o.requires
.PHONY : CMakeFiles/XSModel.dir/requires

CMakeFiles/XSModel.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/XSModel.dir/cmake_clean.cmake
.PHONY : CMakeFiles/XSModel.dir/clean

CMakeFiles/XSModel.dir/depend:
	cd /home/Nickie/JLab/HallA/G4MC/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/Nickie/JLab/HallA/G4MC /home/Nickie/JLab/HallA/G4MC /home/Nickie/JLab/HallA/G4MC/cmake /home/Nickie/JLab/HallA/G4MC/cmake /home/Nickie/JLab/HallA/G4MC/cmake/CMakeFiles/XSModel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/XSModel.dir/depend

