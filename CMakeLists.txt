cmake_minimum_required(VERSION 3.1)
project(SofaShells VERSION 1.0)

find_package(SofaFramework REQUIRED)
find_package(SofaBase REQUIRED)
find_package(SofaGeneral REQUIRED)
find_package(SofaCommon REQUIRED)
find_package(SofaGui QUIET)
find_package(SofaOpenglVisual QUIET)


set(HEADER_FILES

		initPluginShells.h 
		controller/MeshChangedEvent.h 
		controller/MeshInterpolator.h 
		controller/MeshInterpolator.inl 
		controller/TriangleSwitchExample.h 
		controller/TriangleSwitchExample.inl 
		engine/JoinMeshPoints.h 
		engine/JoinMeshPoints.inl 
		engine/FindClosePoints.h 
		engine/FindClosePoints.inl 
		forcefield/BezierTriangularBendingFEMForceField.h 
		forcefield/BezierTriangularBendingFEMForceField.inl 
		forcefield/CstFEMForceField.h 
		forcefield/CstFEMForceField.inl 
		forcefield/TriangularBendingFEMForceField.h 
		forcefield/TriangularBendingFEMForceField.inl 
		forcefield/TriangularShellForceField.h 
		forcefield/TriangularShellForceField.inl 
		mapping/BendingPlateMechanicalMapping.h 
		mapping/BendingPlateMechanicalMapping.inl 
		mapping/BezierTriangleMechanicalMapping.h 
		mapping/BezierTriangleMechanicalMapping.inl 
		shells2/fem/BezierShellInterpolation.h 
		shells2/fem/BezierShellInterpolation.inl 
		shells2/fem/BezierShellInterpolationM.h 
		shells2/fem/BezierShellInterpolationM.inl 
		shells2/forcefield/BezierShellForceField.h 
		shells2/forcefield/BezierShellForceField.inl 
		shells2/mapping/BezierShellMechanicalMapping.h 
		shells2/mapping/BezierShellMechanicalMapping.inl 
		misc/PointProjection.h 
		misc/PointProjection.inl 
	)

set(SOURCE_FILES

		initPluginShells.cpp 
		controller/MeshChangedEvent.cpp
		controller/MeshInterpolator.cpp
		controller/TriangleSwitchExample.cpp 
		engine/JoinMeshPoints.cpp 
		engine/FindClosePoints.cpp 
		forcefield/BezierTriangularBendingFEMForceField.cpp 
		forcefield/CstFEMForceField.cpp 
		forcefield/TriangularBendingFEMForceField.cpp 
		forcefield/TriangularShellForceField.cpp 
		mapping/BendingPlateMechanicalMapping.cpp 
		mapping/BezierTriangleMechanicalMapping.cpp 
		shells2/fem/BezierShellInterpolation.cpp 
		shells2/fem/BezierShellInterpolationM.cpp 
		shells2/forcefield/BezierShellForceField.cpp 
		shells2/mapping/BezierShellMechanicalMapping.cpp 
		misc/PointProjection.cpp

	)
	
	
set(README_FILE "SofaShells.txt" )


option(SOFA-PLUGIN_SHELLS_ADAPTIVITY "Enables shells adaptivity" OFF)


if(SOFA-PLUGIN_SHELLS_ADAPTIVITY)
	set(COMPILER_DEFINE "SOFA_BUILD_SHELLS_ADAPTIVITY")
	
	list(APPEND HEADER_FILES
	
		controller/AdaptiveCuttingController.h 
		controller/AdaptiveCuttingController.inl 
		controller/Test2DAdapter.h 
		controller/Test2DAdapter.inl
		misc/Optimize2DSurface.h 
		misc/Optimize2DSurface.inl 
		misc/SurfaceParametrization.h 
		misc/SurfaceParametrization.inl 
	)


   	list(APPEND SOURCE_FILES
	
		controller/AdaptiveCuttingController.cpp 
		controller/Test2DAdapter.cpp 
		misc/Optimize2DSurface.cpp 
		misc/SurfaceParametrization.cpp
	)


	if(SofaGui_FOUND AND SofaOpenglVisual_FOUND)

		list(APPEND HEADER_FILES
			cutting/AdaptiveCutting.h
			)

		list(APPEND SOURCE_FILES
			cutting/AdaptiveCutting.cpp
			)

	endif()

endif()



add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_SOFASHELLS")

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

target_link_libraries(${PROJECT_NAME} SofaComponentBase SofaMeshCollision SofaGraphComponent SofaCore SofaHelper SofaBaseTopology SofaUserInteraction)
if(SofaGui_FOUND)
	target_link_libraries(${PROJECT_NAME} SofaGuiMain SofaOpenglVisual)
endif()


## Install rules for the library and headers; CMake package configurations files
sofa_generate_package(
    NAME SofaShells VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME}
    INCLUDE_INSTALL_DIR SofaShells
    EXAMPLE_INSTALL_DIR examples
    RELOCATABLE "plugins"
    )

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
SET_PROPERTY(TARGET ${PROJECT_NAME} PROPERTY FOLDER "plugins")