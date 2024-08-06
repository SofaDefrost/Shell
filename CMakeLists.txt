cmake_minimum_required(VERSION 3.12)
project(Shell VERSION 1.0 LANGUAGES CXX)

# Find and load CMake configuration of packages containing this plugin's dependencies
find_package(Sofa.Component.Controller REQUIRED)
sofa_find_package(Sofa.Component.StateContainer REQUIRED)


# List all files
set(SHELL_SRC_DIR src/Shell)
set(HEADER_FILES
    ${SHELL_SRC_DIR}/config.h.in
    ${SHELL_SRC_DIR}/controller/MeshChangedEvent.h
    ${SHELL_SRC_DIR}/controller/MeshInterpolator.h
    ${SHELL_SRC_DIR}/controller/MeshInterpolator.inl
    ${SHELL_SRC_DIR}/engine/JoinMeshPoints.h
    ${SHELL_SRC_DIR}/engine/JoinMeshPoints.inl
    ${SHELL_SRC_DIR}/forcefield/TriangularBendingFEMForceField.h
    ${SHELL_SRC_DIR}/forcefield/TriangularBendingFEMForceField.inl
)
set(SOURCE_FILES
    ${SHELL_SRC_DIR}/initShell.cpp
    ${SHELL_SRC_DIR}/controller/MeshChangedEvent.cpp
    ${SHELL_SRC_DIR}/controller/MeshInterpolator.cpp
    ${SHELL_SRC_DIR}/engine/JoinMeshPoints.cpp
    ${SHELL_SRC_DIR}/forcefield/TriangularBendingFEMForceField.cpp
)
set(README_FILES
    README.md
)

# Create the plugin library.
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

# Link the plugin library to its dependency(ies).
target_link_libraries(${PROJECT_NAME} Sofa.Component.Controller Sofa.Component.StateContainer)

# Create package Config, Version & Target files.
# Deploy the headers, resources, scenes & examples.
# Set the plugin 'relocatable' if built within SOFA.
# --> see SofaMacros.cmake
sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )
