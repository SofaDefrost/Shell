
######  GENERAL PLUGIN CONFIGURATION, you shouldn't have to modify it

SOFA_DIR=../../..
TEMPLATE = lib

include($${SOFA_DIR}/sofa.cfg)

DESTDIR = $$SOFA_DIR/lib/sofa-plugins

#set configuration to dynamic library
CONFIG += $$CONFIGLIBRARIES
CONFIG -= staticlib
CONFIG += dll

###### SPECIFIC PLUGIN CONFIGURATION, you should modify it to configure your plugin

TARGET = PluginShells$$LIBSUFFIX
DEFINES += SOFA_BUILD_SHELLS

LIBS += $$SOFA_LIBS
LIBS += $$SOFA_EXT_LIBS
INCLUDEPATH += $$SOFA_DIR/extlibs

SOURCES = 	forcefield/TriangularBendingFEMForceField.cpp \
			forcefield/TriangularShellForceField.cpp \
			forcefield/BezierTriangularBendingFEMForceField.cpp \
			mapping/BendingPlateMechanicalMapping.cpp \
			mapping/BezierTriangleMechanicalMapping.cpp \
			initPluginShells.cpp

HEADERS = 	forcefield/TriangularBendingFEMForceField.h \
			forcefield/TriangularBendingFEMForceField.inl \
			forcefield/TriangularShellForceField.h \
			forcefield/TriangularShellForceField.inl \
			forcefield/BezierTriangularBendingFEMForceField.h \
			forcefield/BezierTriangularBendingFEMForceField.inl \
			mapping/BendingPlateMechanicalMapping.h \
			mapping/BendingPlateMechanicalMapping.inl \
			mapping/BezierTriangleMechanicalMapping.h \
			mapping/BezierTriangleMechanicalMapping.inl
		  
README_FILE = PluginShells.txt

unix : QMAKE_POST_LINK = cp $$README_FILE $$DESTDIR 
win32 : QMAKE_POST_LINK = copy \"$$README_FILE\" \"$$SOFA_DIR/lib/sofa-plugins\"


