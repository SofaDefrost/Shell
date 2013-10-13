include(${SOFA_CMAKE_DIR}/preProject.cmake)

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
	
	
	
set(README_FILES


	)

option(SOFA-PLUGIN_SHELLS_ADAPTIVITY "Enables shells adaptivity" OFF)

if(SOFA-PLUGIN_SHELLS_ADAPTIVITY)
	set(COMPILER_DEFINE "SOFA_BUILD_SHELLS_ADAPTIVITY")
	
	list(APPEND HEADER_FILES
	
		controller/AdaptiveCuttingController.h 
		controller/AdaptiveCuttingController.inl 
		controller/Test2DAdapter.h 
		controller/Test2DAdapter.inl 
		cutting/AdaptiveCutting.h 
		misc/Optimize2DSurface.h 
		misc/Optimize2DSurface.inl 
		misc/SurfaceParametrization.h 
		misc/SurfaceParametrization.inl 
	)

   	list(APPEND SOURCE_FILES
	
		controller/AdaptiveCuttingController.cpp 
		controller/Test2DAdapter.cpp 
		cutting/AdaptiveCutting.cpp 
		misc/Optimize2DSurface.cpp 
		misc/SurfaceParametrization.cpp
	)	
endif()


# if(LINUX OR APPLE)
# 
# 	include_directories("/usr/include/suitesparse/")
# 
# endif()

include_directories("${SOFA_EXTLIBS_DIR}")
include_directories("${SOFA_APPLICATIONS_PLUGINS_DIR}")

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

AddCompilerDefinitions("SOFA_BUILD_SHELLS")
AddLinkerDependencies(SofaComponentBase SofaHelper SofaBaseTopology SofaUserInteraction)


include(${SOFA_CMAKE_DIR}/postProject.cmake)
