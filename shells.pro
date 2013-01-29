load(sofa/pre)
defineAsPlugin(shells)

TEMPLATE = lib
TARGET = shells


contains (DEFINES, SOFA_QT4) {
        CONFIG += qt
        QT += opengl qt3support xml
}
else {
        CONFIG += qt
        QT += opengl
}


DEFINES += SOFA_BUILD_SHELLS

SOURCES =   initPluginShells.cpp \
            controller/MeshInterpolator.cpp \
            controller/TriangleSwitchExample.cpp \
            engine/JoinMeshPoints.cpp \
            engine/FindClosePoints.cpp \
            forcefield/BezierTriangularBendingFEMForceField.cpp \
            forcefield/TriangularBendingFEMForceField.cpp \
            forcefield/TriangularShellForceField.cpp \
            mapping/BendingPlateMechanicalMapping.cpp \
            mapping/BezierTriangleMechanicalMapping.cpp \
            shells2/fem/BezierShellInterpolation.cpp \
            shells2/forcefield/BezierShellForceField.cpp \
            shells2/mapping/BezierShellMechanicalMapping.cpp \
            misc/PointProjection.cpp \


HEADERS =   initPluginShells.h \
            controller/MeshChangedEvent.h \
            controller/MeshInterpolator.h \
            controller/MeshInterpolator.inl \
            controller/TriangleSwitchExample.h \
            controller/TriangleSwitchExample.inl \
            engine/JoinMeshPoints.h \
            engine/JoinMeshPoints.inl \
            engine/FindClosePoints.h \
            engine/FindClosePoints.inl \
            forcefield/BezierTriangularBendingFEMForceField.h \
            forcefield/BezierTriangularBendingFEMForceField.inl \
            forcefield/TriangularBendingFEMForceField.h \
            forcefield/TriangularBendingFEMForceField.inl \
            forcefield/TriangularShellForceField.h \
            forcefield/TriangularShellForceField.inl \
            mapping/BendingPlateMechanicalMapping.h \
            mapping/BendingPlateMechanicalMapping.inl \
            mapping/BezierTriangleMechanicalMapping.h \
            mapping/BezierTriangleMechanicalMapping.inl \
            shells2/fem/BezierShellInterpolation.h \
            shells2/fem/BezierShellInterpolation.inl \
            shells2/forcefield/BezierShellForceField.h \
            shells2/forcefield/BezierShellForceField.inl \
            shells2/mapping/BezierShellMechanicalMapping.h \
            shells2/mapping/BezierShellMechanicalMapping.inl \
            misc/PointProjection.h \
            misc/PointProjection.inl \


#
# Mesh optimization and adaptivity
#
contains (DEFINES, SOFA_HAVE_SHELL_ADAPTIVITY) {

DEFINES += SOFA_BUILD_SHELLS_ADAPTIVITY

SOURCES += \
            controller/AdaptiveCuttingController.cpp \
            controller/Test2DAdapter.cpp \
            cutting/AdaptiveCutting.cpp \

HEADERS += \
            controller/AdaptiveCuttingController.h \
            controller/AdaptiveCuttingController.inl \
            controller/Test2DAdapter.h \
            controller/Test2DAdapter.inl \
            cutting/AdaptiveCutting.h \
}

# CUDA classes are broken right now.
#contains (DEFINES, SOFA_GPU_CUDA) {
#
#HEADERS += \
#                controller/CudaTest2DAdapter.h \
#                controller/CudaTest2DAdapter.inl
#SOURCES +=      controller/CudaTest2DAdapter.cpp
#CUDA_SOURCES += controller/CudaTest2DAdapter.cu
#
#}

README_FILE = shells.txt

unix : QMAKE_POST_LINK = cp $$SRC_DIR/$$README_FILE $$LIB_DESTDIR
win32 : QMAKE_POST_LINK = copy \"$$toWindowsPath($$SRC_DIR/$$README_FILE)\" \"$$LIB_DESTDIR\"

load(sofa/post)
