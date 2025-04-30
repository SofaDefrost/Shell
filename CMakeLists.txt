cmake_minimum_required(VERSION 3.12)
project(Shell VERSION 1.0 LANGUAGES CXX)

# Find and load CMake configuration of packages containing this plugin's dependencies
find_package(Sofa.Config REQUIRED)
sofa_find_package(Sofa.GL REQUIRED)
sofa_find_package(Sofa.Type REQUIRED)
sofa_find_package(Sofa.Geometry REQUIRED)
sofa_find_package(Sofa.Simulation.Core REQUIRED)
sofa_find_package(Sofa.Component.Controller REQUIRED)
sofa_find_package(Sofa.Component.StateContainer REQUIRED)
sofa_find_package(Sofa.Component.Mapping.NonLinear REQUIRED)
sofa_find_package(Sofa.Component.Mapping.Linear REQUIRED)
sofa_find_package(Sofa.Component.Constraint REQUIRED)
sofa_find_package(Sofa.Component.Constraint.Lagrangian REQUIRED)
sofa_find_package(Sofa.Component.Constraint.Lagrangian.Model REQUIRED)
sofa_find_package(Sofa.Component.Constraint.Lagrangian.Solver REQUIRED)
sofa_find_package(Sofa.Component.Topology.Container.Dynamic REQUIRED)
sofa_find_package(Sofa.Component.Collision.Detection.Intersection REQUIRED)

set(README_FILE README.md)

option(SOFA-PLUGIN_SHELLS_ADAPTIVITY "Enables shells adaptivity" OFF)

# List all files
set(SHELL_SRC_DIR src/Shell)
set(HEADER_FILES
    ${SHELL_SRC_DIR}/config.h.in
    ${SHELL_SRC_DIR}/controller/MeshChangedEvent.h
    ${SHELL_SRC_DIR}/controller/MeshInterpolator.h
    ${SHELL_SRC_DIR}/controller/MeshInterpolator.inl
    ${SHELL_SRC_DIR}/controller/TriangleSwitchExample.h
    ${SHELL_SRC_DIR}/controller/TriangleSwitchExample.inl
    ${SHELL_SRC_DIR}/engine/JoinMeshPoints.h
    ${SHELL_SRC_DIR}/engine/JoinMeshPoints.inl
    ${SHELL_SRC_DIR}/engine/FindClosePoints.h
    ${SHELL_SRC_DIR}/engine/FindClosePoints.inl
    ${SHELL_SRC_DIR}/forcefield/BezierTriangularBendingFEMForceField.h
    ${SHELL_SRC_DIR}/forcefield/BezierTriangularBendingFEMForceField.inl
    ${SHELL_SRC_DIR}/forcefield/CstFEMForceField.h
    ${SHELL_SRC_DIR}/forcefield/CstFEMForceField.inl
    ${SHELL_SRC_DIR}/forcefield/TriangularBendingFEMForceField.h
    ${SHELL_SRC_DIR}/forcefield/TriangularBendingFEMForceField.inl
    ${SHELL_SRC_DIR}/forcefield/TriangularShellForceField.h
    ${SHELL_SRC_DIR}/forcefield/TriangularShellForceField.inl
    ${SHELL_SRC_DIR}/mapping/BendingPlateMechanicalMapping.h
    ${SHELL_SRC_DIR}/mapping/BendingPlateMechanicalMapping.inl
    ${SHELL_SRC_DIR}/mapping/BezierTriangleMechanicalMapping.h
    ${SHELL_SRC_DIR}/mapping/BezierTriangleMechanicalMapping.inl
    ${SHELL_SRC_DIR}/misc/PointProjection.h
    ${SHELL_SRC_DIR}/misc/PointProjection.inl
    ${SHELL_SRC_DIR}/shells2/fem/BezierShellInterpolation.h
    ${SHELL_SRC_DIR}/shells2/fem/BezierShellInterpolation.inl
    ${SHELL_SRC_DIR}/shells2/fem/BezierShellInterpolationM.h
    ${SHELL_SRC_DIR}/shells2/fem/BezierShellInterpolationM.inl
    ${SHELL_SRC_DIR}/shells2/forcefield/BezierShellForceField.h
    ${SHELL_SRC_DIR}/shells2/forcefield/BezierShellForceField.inl
    ${SHELL_SRC_DIR}/shells2/mapping/BezierShellMechanicalMapping.h
    ${SHELL_SRC_DIR}/shells2/mapping/BezierShellMechanicalMapping.inl
)

set(SOURCE_FILES
    ${SHELL_SRC_DIR}/initShell.cpp
    ${SHELL_SRC_DIR}/controller/MeshChangedEvent.cpp
    ${SHELL_SRC_DIR}/controller/MeshInterpolator.cpp
    ${SHELL_SRC_DIR}/controller/TriangleSwitchExample.cpp
    ${SHELL_SRC_DIR}/engine/JoinMeshPoints.cpp
    ${SHELL_SRC_DIR}/engine/FindClosePoints.cpp
    ${SHELL_SRC_DIR}/forcefield/BezierTriangularBendingFEMForceField.cpp
    ${SHELL_SRC_DIR}/forcefield/CstFEMForceField.cpp
    ${SHELL_SRC_DIR}/forcefield/TriangularBendingFEMForceField.cpp
    ${SHELL_SRC_DIR}/forcefield/TriangularShellForceField.cpp
    ${SHELL_SRC_DIR}/mapping/BendingPlateMechanicalMapping.cpp
    ${SHELL_SRC_DIR}/mapping/BezierTriangleMechanicalMapping.cpp
    ${SHELL_SRC_DIR}/misc/PointProjection.cpp
    ${SHELL_SRC_DIR}/shells2/fem/BezierShellInterpolation.cpp
    ${SHELL_SRC_DIR}/shells2/fem/BezierShellInterpolationM.cpp
    ${SHELL_SRC_DIR}/shells2/forcefield/BezierShellForceField.cpp
    ${SHELL_SRC_DIR}/shells2/mapping/BezierShellMechanicalMapping.cpp
)

if(SOFA-PLUGIN_SHELLS_ADAPTIVITY)
    set(COMPILER_DEFINE "SOFA_BUILD_SHELLS_ADAPTIVITY")

    list(APPEND HEADER_FILES
        ${SHELL_SRC_DIR}/controller/AdaptiveCuttingController.h
        ${SHELL_SRC_DIR}/controller/AdaptiveCuttingController.inl
        ${SHELL_SRC_DIR}/controller/Test2DAdapter.h
        ${SHELL_SRC_DIR}/controller/Test2DAdapter.inl
        ${SHELL_SRC_DIR}/misc/Optimize2DSurface.h
        ${SHELL_SRC_DIR}/misc/Optimize2DSurface.inl
        ${SHELL_SRC_DIR}/misc/SurfaceParametrization.h
        ${SHELL_SRC_DIR}/misc/SurfaceParametrization.inl
    )

    list(APPEND SOURCE_FILES
        ${SHELL_SRC_DIR}/controller/AdaptiveCuttingController.cpp
        ${SHELL_SRC_DIR}/controller/Test2DAdapter.cpp
        ${SHELL_SRC_DIR}/misc/Optimize2DSurface.cpp
        ${SHELL_SRC_DIR}/misc/SurfaceParametrization.cpp
    )

    if(SofaGui_FOUND AND SofaOpenglVisual_FOUND)
        list(APPEND HEADER_FILES
            ${SHELL_SRC_DIR}/cutting/AdaptiveCutting.h
        )

        list(APPEND SOURCE_FILES
            ${SHELL_SRC_DIR}/cutting/AdaptiveCutting.cpp
        )
    endif()

endif()


# Create the plugin library
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

# Link the plugin library to its dependency(ies).
target_link_libraries(${PROJECT_NAME}
    Sofa.Component.Controller
    Sofa.Component.Topology.Container.Dynamic
    Sofa.Component.StateContainer
    Sofa.GL
    Sofa.Component.Collision.Detection.Intersection
    Sofa.Component.Mapping.NonLinear
    Sofa.Component.Mapping.Linear
    Sofa.Component.Constraint.Lagrangian.Model
    Sofa.Component.Constraint.Lagrangian.Solver
    Sofa.Component.Constraint.Lagrangian
    Sofa.Component.Constraint.Projective
    Sofa.Component.Mass
    Sofa.Component.SolidMechanics.Spring
    Sofa.Component.MechanicalLoad
    Sofa.Component.LinearSolver.Iterative
    Sofa.Component.ODESolver.Backward
    Sofa.Component.Engine.Select
)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )
