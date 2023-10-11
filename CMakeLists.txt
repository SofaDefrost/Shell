cmake_minimum_required(VERSION 3.12)
project(SofaShells VERSION 1.0 LANGUAGES CXX)


find_package(Sofa.Config REQUIRED)
sofa_find_package(Sofa.Component.Constraint.Lagrangian.Model REQUIRED)
sofa_find_package(Sofa.Component.StateContainer REQUIRED)
sofa_find_package(Sofa.Component.Mapping.NonLinear REQUIRED)
sofa_find_package(Sofa.GL REQUIRED)
sofa_find_package(Sofa.Component.Topology.Container.Dynamic REQUIRED)

set(README_FILE "SofaShells.txt" )
set(SOFASHELLS_SRC_DIR src/SofaShells)


option(SOFA-PLUGIN_SHELLS_ADAPTIVITY "Enables shells adaptivity" OFF)


set(HEADER_FILES

		${SOFASHELLS_SRC_DIR}/config.h.in
                ${SOFASHELLS_SRC_DIR}/controller/MeshChangedEvent.h
                ${SOFASHELLS_SRC_DIR}/controller/MeshInterpolator.h
                ${SOFASHELLS_SRC_DIR}/controller/MeshInterpolator.inl
                ${SOFASHELLS_SRC_DIR}/controller/TriangleSwitchExample.h
                ${SOFASHELLS_SRC_DIR}/controller/TriangleSwitchExample.inl
                ${SOFASHELLS_SRC_DIR}/engine/JoinMeshPoints.h
                ${SOFASHELLS_SRC_DIR}/engine/JoinMeshPoints.inl
                ${SOFASHELLS_SRC_DIR}/engine/FindClosePoints.h
                ${SOFASHELLS_SRC_DIR}/engine/FindClosePoints.inl
                ${SOFASHELLS_SRC_DIR}/forcefield/BezierTriangularBendingFEMForceField.h
                ${SOFASHELLS_SRC_DIR}/forcefield/BezierTriangularBendingFEMForceField.inl
                ${SOFASHELLS_SRC_DIR}/forcefield/CstFEMForceField.h
                ${SOFASHELLS_SRC_DIR}/forcefield/CstFEMForceField.inl
                ${SOFASHELLS_SRC_DIR}/forcefield/TriangularBendingFEMForceField.h
                ${SOFASHELLS_SRC_DIR}/forcefield/TriangularBendingFEMForceField.inl
                ${SOFASHELLS_SRC_DIR}/forcefield/TriangularShellForceField.h
                ${SOFASHELLS_SRC_DIR}/forcefield/TriangularShellForceField.inl
                ${SOFASHELLS_SRC_DIR}/mapping/BendingPlateMechanicalMapping.h
                ${SOFASHELLS_SRC_DIR}/mapping/BendingPlateMechanicalMapping.inl
                ${SOFASHELLS_SRC_DIR}/mapping/BezierTriangleMechanicalMapping.h
                ${SOFASHELLS_SRC_DIR}/mapping/BezierTriangleMechanicalMapping.inl
                ${SOFASHELLS_SRC_DIR}/misc/PointProjection.h
                ${SOFASHELLS_SRC_DIR}/misc/PointProjection.inl
                ${SOFASHELLS_SRC_DIR}/shells2/fem/BezierShellInterpolation.h
                ${SOFASHELLS_SRC_DIR}/shells2/fem/BezierShellInterpolation.inl
                ${SOFASHELLS_SRC_DIR}/shells2/fem/BezierShellInterpolationM.h
                ${SOFASHELLS_SRC_DIR}/shells2/fem/BezierShellInterpolationM.inl
                ${SOFASHELLS_SRC_DIR}/shells2/forcefield/BezierShellForceField.h
                ${SOFASHELLS_SRC_DIR}/shells2/forcefield/BezierShellForceField.inl
                ${SOFASHELLS_SRC_DIR}/shells2/mapping/BezierShellMechanicalMapping.h
                ${SOFASHELLS_SRC_DIR}/shells2/mapping/BezierShellMechanicalMapping.inl
	)

set(SOURCE_FILES

		${SOFASHELLS_SRC_DIR}/initPluginShells.cpp 
		${SOFASHELLS_SRC_DIR}/controller/MeshChangedEvent.cpp
                ${SOFASHELLS_SRC_DIR}/controller/MeshInterpolator.cpp
                ${SOFASHELLS_SRC_DIR}/controller/TriangleSwitchExample.cpp
                ${SOFASHELLS_SRC_DIR}/engine/JoinMeshPoints.cpp
                ${SOFASHELLS_SRC_DIR}/engine/FindClosePoints.cpp
                ${SOFASHELLS_SRC_DIR}/forcefield/BezierTriangularBendingFEMForceField.cpp
                ${SOFASHELLS_SRC_DIR}/forcefield/CstFEMForceField.cpp
                ${SOFASHELLS_SRC_DIR}/forcefield/TriangularBendingFEMForceField.cpp
                ${SOFASHELLS_SRC_DIR}/forcefield/TriangularShellForceField.cpp
                ${SOFASHELLS_SRC_DIR}/mapping/BendingPlateMechanicalMapping.cpp
                ${SOFASHELLS_SRC_DIR}/mapping/BezierTriangleMechanicalMapping.cpp
                ${SOFASHELLS_SRC_DIR}/misc/PointProjection.cpp
                ${SOFASHELLS_SRC_DIR}/shells2/fem/BezierShellInterpolation.cpp
                ${SOFASHELLS_SRC_DIR}/shells2/fem/BezierShellInterpolationM.cpp
                ${SOFASHELLS_SRC_DIR}/shells2/forcefield/BezierShellForceField.cpp
                ${SOFASHELLS_SRC_DIR}/shells2/mapping/BezierShellMechanicalMapping.cpp
	)
	
	




if(SOFA-PLUGIN_SHELLS_ADAPTIVITY)
	set(COMPILER_DEFINE "SOFA_BUILD_SHELLS_ADAPTIVITY")
	
	list(APPEND HEADER_FILES
	
		${SOFASHELLS_SRC_DIR}/controller/AdaptiveCuttingController.h 
		${SOFASHELLS_SRC_DIR}/controller/AdaptiveCuttingController.inl 
		${SOFASHELLS_SRC_DIR}/controller/Test2DAdapter.h 
		${SOFASHELLS_SRC_DIR}/controller/Test2DAdapter.inl
		${SOFASHELLS_SRC_DIR}/misc/Optimize2DSurface.h 
		${SOFASHELLS_SRC_DIR}/misc/Optimize2DSurface.inl 
		${SOFASHELLS_SRC_DIR}/misc/SurfaceParametrization.h 
		${SOFASHELLS_SRC_DIR}/misc/SurfaceParametrization.inl 
	)


   	list(APPEND SOURCE_FILES
	
		${SOFASHELLS_SRC_DIR}/controller/AdaptiveCuttingController.cpp 
		${SOFASHELLS_SRC_DIR}/controller/Test2DAdapter.cpp 
		${SOFASHELLS_SRC_DIR}/misc/Optimize2DSurface.cpp 
		${SOFASHELLS_SRC_DIR}/misc/SurfaceParametrization.cpp
	)


	if(SofaGui_FOUND AND SofaOpenglVisual_FOUND)

		list(APPEND HEADER_FILES
			${SOFASHELLS_SRC_DIR}/cutting/AdaptiveCutting.h
			)

		list(APPEND SOURCE_FILES
			${SOFASHELLS_SRC_DIR}/cutting/AdaptiveCutting.cpp
			)

	endif()

endif()


# Create the plugin library
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

# Link the plugin library to its dependency
target_link_libraries(${PROJECT_NAME} SofaCore SofaHelper SofaBaseCollision SofaBaseTopology SofaBaseLinearSolver SofaBaseMechanics SofaUserInteraction SofaMeshCollision SofaGraphComponent SofaBoundaryCondition Sofa.GL)


if(SofaGui_FOUND AND SofaOpenglVisual_FOUND)
	target_link_libraries(${PROJECT_NAME} SofaGui SofaOpenglVisual)
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
