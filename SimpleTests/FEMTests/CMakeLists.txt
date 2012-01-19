
#FIND_PACKAGE ( Eigen REQUIRED )
#FIND_PACKAGE ( LIS REQUIRED )

IF(TEST_Matrix)
IF(NOT MSVC)
	## pthread ##
	SET ( CMAKE_THREAD_PREFER_PTHREAD On )
	FIND_PACKAGE( Threads )

	# Find blas
	FIND_PACKAGE ( BLAS REQUIRED )
	# Find lapack
	FIND_PACKAGE ( LAPACK REQUIRED )
        IF ( CMAKE_USE_PTHREADS_INIT )
        	SET (HAVE_PTHREADS TRUE)
        	MESSAGE (STATUS "pthread library found." )
	ENDIF (CMAKE_USE_PTHREADS_INIT )
ENDIF(NOT MSVC)
ENDIF(TEST_Matrix)


INCLUDE_DIRECTORIES(
        .
	../../Base/
	../../MathLib/
	../../GeoLib/
	../../MeshLib/
	../../FemLib/
)

# Create the executable
ADD_EXECUTABLE( FEMTests
        FEMTests.cpp
        ${SOURCES}
        ${HEADERS}
)

#IF (WIN32)
#        TARGET_LINK_LIBRARIES(MatMult Winmm.lib)
#ENDIF (WIN32)

IF (LIS_FOUND)
	TARGET_LINK_LIBRARIES( FEMTests ${LIS_LIBRARIES} )
ENDIF (LIS_FOUND)
IF(NOT MSVC)
	TARGET_LINK_LIBRARIES ( FEMTests ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
ENDIF(NOT MSVC)


TARGET_LINK_LIBRARIES ( FEMTests
	Base
	MathLib
	GeoLib
#	MeshLib
#	FemLib
)

###WW
IF(USE_PETSC)
        target_link_libraries(FEMTests  ${PETSC_LIBRARIES})
ENDIF(USE_PETSC)
