# CMakeLists.txt for Shell.Adaptivity plugin
cmake_minimum_required(VERSION 3.12)
project(Shell.Adaptivity VERSION 1.0 LANGUAGES CXX)

find_package(Sofa.Config REQUIRED)
sofa_find_package(Shell REQUIRED)

set(HEADER_FILES
    src/Shell/Adaptivity/config.h.in
    src/Shell/Adaptivity/init.h

    src/Shell/Adaptivity/controller/AdaptiveCuttingController.h
    src/Shell/Adaptivity/controller/AdaptiveCuttingController.inl
    src/Shell/Adaptivity/controller/Test2DAdapter.h
    src/Shell/Adaptivity/controller/Test2DAdapter.inl

    src/Shell/Adaptivity/cutting/AdaptiveCutting.h

    src/Shell/Adaptivity/misc/Optimize2DSurface.h
    src/Shell/Adaptivity/misc/Optimize2DSurface.inl
    src/Shell/Adaptivity/misc/SurfaceParametrization.h
    src/Shell/Adaptivity/misc/SurfaceParametrization.inl
)

set(SOURCE_FILES
    src/Shell/Adaptivity/init.cpp

    src/Shell/Adaptivity/controller/AdaptiveCuttingController.cpp
    src/Shell/Adaptivity/controller/Test2DAdapter.cpp

    src/Shell/Adaptivity/cutting/AdaptiveCutting.h

    src/Shell/Adaptivity/misc/Optimize2DSurface.cpp
    src/Shell/Adaptivity/misc/SurfaceParametrization.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Shell)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )
