def SetupEigenstates1D(prop, potentialIndices=[0]):
	"""
	Finds the eigenvalues and eigenvectors of the first potential
	of prop. From the default config file, this is the field free
	1D Hydrogen system.

	eigenvalues is a list of 1-d eigenvalue arrays. Each array corresponding
	to a
	"""
	if not pyprop.IsSingleProc():
		raise Exception("Works only on a single processor")

	eigenValues = []
	eigenVectors = []

	eigenvectorScaling = 1

	M = SetupPotentialMatrix(prop, potentialIndices)

	E, V = scipy.linalg.eig(M)

	idx = argsort(real(E))
	E = real(E[idx])
	eigenValues = E

	#Sort and normalize eigenvectors
	sNorm = lambda v: sqrt(abs(sum(conj(v) * dot(S, v))))
	eigenVectors = array([v/sNorm(v) for v in [V[:,idx[i]] for i in range(V.shape[1])]]).transpose()

	#assure correct phase convention (first oscillation should start out real positive)
	for i, curE in enumerate(E):
		bspl.ConstructFunctionFromBSplineExpansion(V[:,i].copy(), phaseGrid, phaseBuffer)
		phase = arctan2(imag(phaseBuffer[1]), real(phaseBuffer[1]))
		eigenVectors[:,i] *= exp(-1.0j * phase)

	return eigenValues, eigenVectors


def SaveEigenstates1D(**args):
	postfix = GetGridPostfix(**args)
	outputFile = "eigenstates/eigenstates_1D_%s.h5" % ("_".join(postfix))

	#Setup problem
	prop = SetupProblem(**args)

	#Setp eigenvalues and eigenvectors
	eigenValues, eigenVectors = SetupEigenstates1D(prop)

	#Save eigenvalues and eigenvectors to file
	if outputFile != None:
		f = tables.openFile(outputFile, "w")
		try:
			#Create main group
			eigGroup = f.createGroup(f.root, "Eig")

			#save config object
			eigGroup._v_attrs.configObject = prop.Config.cfgObj

			#save eigenvalues and eigenstates
			f.createArray(eigGroup, "eigenvalues", eigenValues)
			f.createArray(eigGroup, "eigenvectors", eigenVectors)

		finally:
			f.close()


def SetupPotentialMatrix(prop, whichPotentials):
	matrixSize = prop.psi.GetData().shape[0]
	matrix = zeros((matrixSize, matrixSize), dtype=complex)

	for potNum in whichPotentials:	
		if isinstance(potNum, pyprop.TensorPotential):
			potential = potNum
		else:
			potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential: %s" % (potential.Name, )

		basisPairs = potential.BasisPairs[0]

		for i, (x,xp) in enumerate(basisPairs):
			indexLeft = x
			indexRight = xp
			matrix[indexLeft, indexRight] += potential.PotentialData[i]

	return matrix

