"""
OneElectronTools
================

Toolkit for 1D/1-particle related tasks

"""
import sys
import logging
import numpy
import scipy
from scipy import linalg
import tables

from utils import RegisterAll, GetClassLogger
from example import SetupProblem, GetGridPrefix, GetEigenstateFilename
import pyprop


@RegisterAll
def GetHamiltonMatrix1D(prop):
	"""
	Calculate matrix representation of Hamiltonian given by the
	pyprop.Problem-object 'prop'. Only works for 1D problems.

	"""
	#Only for 1D problems
	assert numpy.rank(prop.psi.GetData()) == 1

	size = prop.psi.GetData().size
	matrix = numpy.zeros((size, size), dtype=complex)
	tempPsi = prop.psi.CopyDeep()

	for i in range(size):
		prop.psi.Clear()
		prop.psi.GetData()[i] = 1

		tempPsi.Clear()
		prop.MultiplyHamiltonian(prop.psi, tempPsi)
		
		matrix[:, i] = tempPsi.GetData()[:]
		
	return matrix

@RegisterAll
def SetupEigenstates1D(prop):
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

	#get logger
	#logger = logging.getLogger("%s" % (sys._getframe().f_code.co_name))
	logger = logging.getLogger("%s" % __name__)

	#Get matrix representation of hamiltonian
	logger.info("Setting up hamiltonian matrix representation...")
	M = GetHamiltonMatrix1D(prop)

	#compute eigenvalues
	logger.info("Computing eigenvalues...")
	E, V = scipy.linalg.eigh(M.real)

	#Sort and normalize eigenvectors
	logger.info("Sort and normalize...")
	idx = numpy.argsort(E)
	E = E[idx]
	eigenValues = E
	eigenVectors = numpy.array([v/linalg.norm(v) 
		for v in [V[:,idx[i]] for i in range(V.shape[1])]]).transpose()

	#assure correct phase convention (first oscillation should start out real positive)
	logger.info("Applying phase convention...")
	for i, curE in enumerate(E):
		phaseBuffer = V[:,i]
		phase = numpy.arctan2(numpy.imag(phaseBuffer[1]), numpy.real(phaseBuffer[1]))
		eigenVectors[:,i] *= numpy.exp(-1.0j * phase)

	return eigenValues, eigenVectors


@RegisterAll
def SaveEigenstates1D(**args):
	#prefix = GetGridPrefix(**args)
	#outputFile = "%s/eigenstates_1D_%s.h5" % (EIGENSTATE_PATH, prefix)
	outputFile = GetEigenstateFilename(**args)

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

