pluginList=['SofaEngine', 'SofaLoader', 'SofaSimpleFem','SofaExporter', 'SofaPython3', 'SofaShells']

def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=pluginList)
    rootNode.gravity = [0, -98.1, 0]
    rootNode.addObject('AttachBodyButtonSetting', stiffness=0.1)
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', maxIterations=1e3, tolerance=1e-3)


    square = rootNode.addChild('Square')
    square.addObject('EulerImplicitSolver')
    square.addObject('SparseLDLSolver')
    square.addObject('GenericConstraintCorrection')
    square.addObject('MeshOBJLoader', filename='mesh/square1.obj')
    square.addObject('MeshTopology', src=square.MeshOBJLoader.getLinkPath(), name="")
    square.addObject('MechanicalObject', template='Rigid3', showObject=1)
    square.addObject('UniformMass', totalMass=0.005)
    square.addObject('BoxROI', box=[0, 0.9, -0.1, 1, 1, 0.1], name="BoxROI", drawBoxes=True)
    square.addObject('FixedConstraint', indices=square.BoxROI.indices.getLinkPath())
    square.addObject('TriangularShellForceField', youngModulus=1.7e3, poissonRatio=0.3, thickness=0.05, measure="Von Mises stress")
    visu = square.addChild('Visu')
    visu.addObject('OglModel', src=square.MeshTopology.getLinkPath())
    visu.addObject('IdentityMapping')

