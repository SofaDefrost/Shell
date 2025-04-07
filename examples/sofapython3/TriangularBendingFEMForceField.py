def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', name="Shell")
    rootNode.gravity = [0, -98.1, 0]
    rootNode.addObject('AttachBodyButtonSetting', stiffness=0.1)
    rootNode.addObject('DefaultAnimationLoop')

    square = rootNode.addChild('Square')
    square.addObject('EulerImplicitSolver')
    square.addObject('SparseLDLSolver')
    square.addObject('MeshOBJLoader', filename='mesh/square1.obj')
    square.addObject('MeshTopology', src=square.MeshOBJLoader.getLinkPath())
    square.addObject('MechanicalObject', template='Rigid3')
    square.addObject('UniformMass', totalMass=0.005)
    square.addObject('BoxROI', box=[0, 0.9, -0.1, 1, 1, 0.1], drawBoxes=True)
    square.addObject('FixedConstraint', indices=square.BoxROI.indices.getLinkPath())
    square.addObject('TriangularBendingFEMForceField', youngModulus=1.7e3,
                     poissonRatio=0.3, thickness=0.01)

    visu = square.addChild('Visu')
    visu.addObject('OglModel', src=square.MeshTopology.getLinkPath())
    visu.addObject('IdentityMapping')

