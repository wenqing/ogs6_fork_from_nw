#####################
### OGS-6 Project ###
#####################

# Specify minimum CMake version
cmake_minimum_required(VERSION 2.6)

# Project name
project( OGS-6 )

### CMake includes ###
INCLUDE(scripts/cmake/CMakeSetup.cmake)
INCLUDE(scripts/cmake/CompilerSetup.cmake)
INCLUDE(scripts/cmake/Find.cmake)

###############
### Options ###
###############
###
###WW
OPTION(USE_PETSC "Use PETSc libs for parallel computing" OFF)
OPTION(USE_LIS "Use LIS lib" OFF)
OPTION(USE_EIGEN "Use EIGEN for local assembly" OFF)
OPTION(USE_OPENMP "Use OpenMP" OFF)
#
OPTION(TEST_Solver "Solver test" OFF)
OPTION(TEST_Matrix "Matrix test" OFF)
OPTION(TEST_FEM "FEM test" ON)


# Profiling
IF((CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC) AND GPROF_PATH)
	OPTION(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
ENDIF() # GCC AND GPROF_PATH


IF(USE_OPENMP)
  ADD_DEFINITIONS( _OPENMP)
ENDIF(USE_OPENMP)

IF(USE_EIGEN)
  FIND_PACKAGE ( Eigen REQUIRED )
  ADD_DEFINITIONS(-DUSE_PETSC)
ENDIF(USE_EIGEN)

IF(USE_LIS)
  FIND_PACKAGE ( LIS REQUIRED )
   ADD_DEFINITIONS(-DUSE_LIS)
ENDIF(USE_LIS)

#### For PETSc
###WW
IF(USE_PETSC)
MESSAGE (STATUS  "Configuring for PETSc" )

	ADD_DEFINITIONS(-DUSE_PETSC
	  -Wunused-parameter
        )


####	SET(FIND_PETSC ${PROJECT_SOURCE_DIR})
	SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake )


	FIND_PACKAGE(PETSc REQUIRED)


       include_directories(
              ${PETSC_INCLUDES} 
        )

	FIND_PACKAGE(MPI)
	IF(MPI_FOUND)		
	  SET(CMAKE_C_COMPILER ${MPI_COMPILER})
	  SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
	ELSE(MPI_FOUND)
	  MESSAGE (FATAL_ERROR "Aborting: MPI implementation is not found!")
	ENDIF(MPI_FOUND)			
ENDIF(USE_PETSC)



# Set build directories
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

# Add subdirectories with the projects
ADD_SUBDIRECTORY( Base )
ADD_SUBDIRECTORY( GeoLib )
ADD_SUBDIRECTORY( MathLib )
ADD_SUBDIRECTORY( MeshLib )
ADD_SUBDIRECTORY( FemLib )

IF(TEST_Solver)
ADD_SUBDIRECTORY( SimpleTests/MatrixTests )
ENDIF(TEST_Solver)

IF(TEST_Matrix)
	IF(NOT MSVC)
	        ADD_SUBDIRECTORY( SimpleTests/SolverTests )
	ENDIF(NOT MSVC)
ENDIF(TEST_Matrix)

IF(TEST_FEM)
	ADD_SUBDIRECTORY( SimpleTests/FEMTests )
ENDIF(TEST_FEM)

