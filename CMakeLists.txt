cmake_minimum_required(VERSION 3.12)
project(SofaShells VERSION 1.0 LANGUAGES CXX)

find_package(SofaFramework REQUIRED)
find_package(SofaBase REQUIRED)
find_package(SofaGeneral REQUIRED)
find_package(SofaCommon REQUIRED)
find_package(SofaGui QUIET)
find_package(SofaOpenglVisual QUIET)

set(README_FILE "SofaShells.txt" )
set(SOFASHELLS_SRC_DIR src/SofaShells)

option(SOFA-PLUGIN_SHELLS_ADAPTIVITY "Enables shells adaptivity" OFF)


set(HEADER_FILES

		${SOFAMJEDFEM_SRC_DIR}/initPluginShells.h 
		${SOFAMJEDFEM_SRC_DIR}/controller/MeshChangedEvent.h 
		${SOFAMJEDFEM_SRC_DIR}/controller/MeshInterpolator.h 
		${SOFAMJEDFEM_SRC_DIR}/controller/MeshInterpolator.inl 
		${SOFAMJEDFEM_SRC_DIR}/controller/TriangleSwitchExample.h 
		${SOFAMJEDFEM_SRC_DIR}/controller/TriangleSwitchExample.inl 
		${SOFAMJEDFEM_SRC_DIR}/engine/JoinMeshPoints.h 
		${SOFAMJEDFEM_SRC_DIR}/engine/JoinMeshPoints.inl 
		${SOFAMJEDFEM_SRC_DIR}/engine/FindClosePoints.h 
		${SOFAMJEDFEM_SRC_DIR}/engine/FindClosePoints.inl 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/BezierTriangularBendingFEMForceField.h 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/BezierTriangularBendingFEMForceField.inl 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/CstFEMForceField.h 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/CstFEMForceField.inl 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/TriangularBendingFEMForceField.h 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/TriangularBendingFEMForceField.inl 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/TriangularShellForceField.h 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/TriangularShellForceField.inl 
		${SOFAMJEDFEM_SRC_DIR}/mapping/BendingPlateMechanicalMapping.h 
		${SOFAMJEDFEM_SRC_DIR}/mapping/BendingPlateMechanicalMapping.inl 
		${SOFAMJEDFEM_SRC_DIR}/mapping/BezierTriangleMechanicalMapping.h 
		${SOFAMJEDFEM_SRC_DIR}/mapping/BezierTriangleMechanicalMapping.inl 
		${SOFAMJEDFEM_SRC_DIR}/shells2/fem/BezierShellInterpolation.h 
		${SOFAMJEDFEM_SRC_DIR}/shells2/fem/BezierShellInterpolation.inl 
		${SOFAMJEDFEM_SRC_DIR}/shells2/fem/BezierShellInterpolationM.h 
		${SOFAMJEDFEM_SRC_DIR}/shells2/fem/BezierShellInterpolationM.inl 
		${SOFAMJEDFEM_SRC_DIR}/shells2/forcefield/BezierShellForceField.h 
		${SOFAMJEDFEM_SRC_DIR}/shells2/forcefield/BezierShellForceField.inl 
		${SOFAMJEDFEM_SRC_DIR}/shells2/mapping/BezierShellMechanicalMapping.h 
		${SOFAMJEDFEM_SRC_DIR}/shells2/mapping/BezierShellMechanicalMapping.inl 
		${SOFAMJEDFEM_SRC_DIR}/misc/PointProjection.h 
		${SOFAMJEDFEM_SRC_DIR}/misc/PointProjection.inl 
	)

set(SOURCE_FILES

		${SOFAMJEDFEM_SRC_DIR}/initPluginShells.cpp 
		${SOFAMJEDFEM_SRC_DIR}/controller/MeshChangedEvent.cpp
		${SOFAMJEDFEM_SRC_DIR}/controller/MeshInterpolator.cpp
		${SOFAMJEDFEM_SRC_DIR}/controller/TriangleSwitchExample.cpp 
		${SOFAMJEDFEM_SRC_DIR}/engine/JoinMeshPoints.cpp 
		${SOFAMJEDFEM_SRC_DIR}/engine/FindClosePoints.cpp 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/BezierTriangularBendingFEMForceField.cpp 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/CstFEMForceField.cpp 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/TriangularBendingFEMForceField.cpp 
		${SOFAMJEDFEM_SRC_DIR}/forcefield/TriangularShellForceField.cpp 
		${SOFAMJEDFEM_SRC_DIR}/mapping/BendingPlateMechanicalMapping.cpp 
		${SOFAMJEDFEM_SRC_DIR}/mapping/BezierTriangleMechanicalMapping.cpp 
		${SOFAMJEDFEM_SRC_DIR}/shells2/fem/BezierShellInterpolation.cpp 
		${SOFAMJEDFEM_SRC_DIR}/shells2/fem/BezierShellInterpolationM.cpp 
		${SOFAMJEDFEM_SRC_DIR}/shells2/forcefield/BezierShellForceField.cpp 
		${SOFAMJEDFEM_SRC_DIR}/shells2/mapping/BezierShellMechanicalMapping.cpp 
		${SOFAMJEDFEM_SRC_DIR}/misc/PointProjection.cpp
	)
	
	




if(SOFA-PLUGIN_SHELLS_ADAPTIVITY)
	set(COMPILER_DEFINE "SOFA_BUILD_SHELLS_ADAPTIVITY")
	
	list(APPEND HEADER_FILES
	
		${SOFAMJEDFEM_SRC_DIR}/controller/AdaptiveCuttingController.h 
		${SOFAMJEDFEM_SRC_DIR}/controller/AdaptiveCuttingController.inl 
		${SOFAMJEDFEM_SRC_DIR}/controller/Test2DAdapter.h 
		${SOFAMJEDFEM_SRC_DIR}/controller/Test2DAdapter.inl
		${SOFAMJEDFEM_SRC_DIR}/misc/Optimize2DSurface.h 
		${SOFAMJEDFEM_SRC_DIR}/misc/Optimize2DSurface.inl 
		${SOFAMJEDFEM_SRC_DIR}/misc/SurfaceParametrization.h 
		${SOFAMJEDFEM_SRC_DIR}/misc/SurfaceParametrization.inl 
	)


   	list(APPEND SOURCE_FILES
	
		${SOFAMJEDFEM_SRC_DIR}/controller/AdaptiveCuttingController.cpp 
		${SOFAMJEDFEM_SRC_DIR}/controller/Test2DAdapter.cpp 
		${SOFAMJEDFEM_SRC_DIR}/misc/Optimize2DSurface.cpp 
		${SOFAMJEDFEM_SRC_DIR}/misc/SurfaceParametrization.cpp
	)


	if(SofaGui_FOUND AND SofaOpenglVisual_FOUND)

		list(APPEND HEADER_FILES
			${SOFAMJEDFEM_SRC_DIR}/cutting/AdaptiveCutting.h
			)

		list(APPEND SOURCE_FILES
			${SOFAMJEDFEM_SRC_DIR}/cutting/AdaptiveCutting.cpp
			)

	endif()

endif()


# Create the plugin library
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

# Link the plugin library to its dependency
target_link_libraries(${PROJECT_NAME} SofaComponentBase SofaMeshCollision SofaGraphComponent SofaCore SofaHelper SofaBaseTopology SofaUserInteraction)

if(SofaGui_FOUND)
	target_link_libraries(${PROJECT_NAME} SofaGuiMain SofaOpenglVisual)
endif()


# Create package Config, Version & Target files
sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )