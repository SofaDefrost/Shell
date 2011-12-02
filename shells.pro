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

SOURCES = 	forcefield/BezierTriangularBendingFEMForceField.cpp \
			forcefield/TriangularBendingFEMForceField.cpp \
			forcefield/TriangularShellForceField.cpp \
			mapping/BendingPlateMechanicalMapping.cpp \
			#mapping/BezierTriangleMechanicalMapping.cpp \
			initPluginShells.cpp

HEADERS =   initPluginShells.h \
			forcefield/BezierTriangularBendingFEMForceField.h \
			forcefield/BezierTriangularBendingFEMForceField.inl \
			forcefield/TriangularBendingFEMForceField.h \
			forcefield/TriangularBendingFEMForceField.inl \
			forcefield/TriangularShellForceField.h \
			forcefield/TriangularShellForceField.inl \
			mapping/BendingPlateMechanicalMapping.h \
			mapping/BendingPlateMechanicalMapping.inl \
			#mapping/BezierTriangleMechanicalMapping.h \
			#mapping/BezierTriangleMechanicalMapping.inl
		  
README_FILE = shells.txt

unix : QMAKE_POST_LINK = cp $$SRC_DIR/$$README_FILE $$LIB_DESTDIR
win32 : QMAKE_POST_LINK = copy \"$$toWindowsPath($$SRC_DIR/$$README_FILE)\" \"$$LIB_DESTDIR\"

load(sofa/post)
