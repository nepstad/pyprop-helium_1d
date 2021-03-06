import os, sys
import pypar
import pyprop
import tables
import numpy
import logging

from helium_1d import MODULE_PATH
from utils import RegisterProjectNamespace, RegisterAll
from laserfunctions import LaserFunction
import copy as objectcopier

#Data file locations. Should be updated by user when the module is imported
EIGENSTATE_PATH = "%s/eigenstates" % MODULE_PATH
GROUNDSTATE_PATH = "%s/groundstate" % MODULE_PATH
SIMULATION_DATA_PATH = "%s/output" % MODULE_PATH

@RegisterProjectNamespace
def SetupConfig(**args):
	configFile = args.get("config", "config.ini")
	#if configfile is a string, load it, otherwise, treat it as
	#a config parser object
	if isinstance(configFile, str):
		conf = pyprop.Load(configFile)
	elif isinstance(configFile, pyprop.Config):
		#Let's make a deep copy here to avoid changing input
		conf = objectcopier.deepcopy(configFile)
	else:
		conf = pyprop.Config(configFile)

	if "imtime" in args:
		if args["imtime"]:
			conf.SetValue("Propagation", "timestep", -1.0j * abs(conf.Propagation.timestep))
			conf.SetValue("Propagation", "renormalization", True)
		else:
			conf.SetValue("Propagation", "timestep", abs(conf.Propagation.timestep))
			conf.SetValue("Propagation", "renormalization", False)

	if "additional_potentials" in args:
		curPot = conf.Propagation.potential_evaluation
		newPot = curPot + args["additional_potentials"]
		conf.SetValue("Propagation", "potential_evaluation", newPot)

	return conf


@RegisterAll
def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	prop.psi.Normalize()

	return prop


@RegisterAll
def GetGridPrefix(**args):
	conf = SetupConfig(**args)
	xmax = conf.Representation.rank0[1]
	count = conf.Representation.rank0[2]
	return "xmax%s_n%s" % (xmax, count)

@RegisterAll
def GetEigenstateFilename(**args):
	prefix = GetGridPrefix(**args)
	filename = "%s/eigenstates_1D_%s.h5" % (EIGENSTATE_PATH, prefix)
	return filename

@RegisterAll
def GetEigenvectorDatasetPath(eigenvectorIndex):
	return "/Eig/%s" % GetEigenvectorDatasetName(eigenvectorIndex) 

@RegisterAll
def GetGroundstateFilename(**args):
	prefix = GetGridPrefix(**args)
	modelPostfix = "_" + args.get("modelPostfix", "std")
	filename = "%s/groundstate_%s%s.h5" % (GROUNDSTATE_PATH, prefix, modelPostfix)
	return filename


@RegisterAll
def FindGroundstate(**args):
	"""
	Use imaginary time propagation to find ground state
	"""
	args["imtime"] = True

	#filename = "%s/groundstate_%s.h5" % (outputDir, GetGridPrefix(**args))
	filename = GetGroundstateFilename(**args)
	dir = os.path.dirname(filename)
	if not os.path.exists(dir) and pyprop.ProcId == 0:
		os.makedirs(dir)

	#Setup problem
	prop = SetupProblem(**args)

	#Propagate in imaginary time to obtain ground state
	for t in prop.Advance(10):
		N = prop.psi.GetNorm()
		E = prop.GetEnergy()
		pyprop.PrintOut("t = %03f, N = %s, E = %s" % (t, N, E))

	#Save groundstate
	prop.psi.Normalize()
	prop.SaveWavefunctionHDF(filename, "/wavefunction")

	return prop


@RegisterAll
def FindEigenvalues(**args):
	"""
	Use Piram to find some eigenvalues.

	"""
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	solver.Solve()
	print solver.GetEigenvalues()

	return solver


@RegisterAll
def SaveEigenvalues(filename, solver):
	""" Saves the output of FindEigenvalues to a hdf5 file.

	"""

	#Get eigenvalues
	prop = solver.BaseProblem
	E = numpy.array(solver.GetEigenvalues())

	#Get eigenvalue error estimates
	errorEstimatesPIRAM = solver.Solver.GetErrorEstimates()
	convergenceEstimatesEig = solver.Solver.GetConvergenceEstimates()

	#remove file if it exists
	try:
		if os.path.exists(filename):
			if pyprop.ProcId == 0:
				os.remove(filename)
	except:
		logging.error("Could not remove %s (%s)" % (filename, sys.exc_info()[1]))

	#Store eigenvalues and eigenvectors
	logging.info("Now storing eigenvectors...")
	for i in range(len(E)):
		solver.SetEigenvector(prop.psi, i)
		prop.SaveWavefunctionHDF(filename, GetEigenvectorDatasetPath(i))

	if pyprop.ProcId == 0:
		pyprop.RemoveExistingDataset(filename, "/Eig/Eigenvalues")
		pyprop.RemoveExistingDataset(filename, "/Eig/ErrorEstimateListPIRAM")
		pyprop.RemoveExistingDataset(filename, "/Eig/ConvergenceEstimateEig")
		h5file = tables.openFile(filename, "r+")
		try:
			myGroup = h5file.getNode("/Eig")
			h5file.createArray(myGroup, "Eigenvalues", E)
			h5file.createArray(myGroup, "ErrorEstimateListPIRAM", errorEstimatesPIRAM)
			h5file.createArray(myGroup, "ConvergenceEstimateEig", convergenceEstimatesEig)

			#Store config
			myGroup._v_attrs.configObject = prop.Config.cfgObj
			
		finally:
			h5file.close()


@RegisterAll
def Propagate(**args):
	#args['imtime'] = False
	#args['additional_potentials'] = ["LaserPotential"]
	#args['config'] = "propagation.ini"
	numOutput = args.get("numOutput", 50)

	#get propagation tasks
	propTasks = args.get("propagationTasks", [])
	
	#Setup problem
	prop = SetupProblem(**args)

	#Get initial state
	#filename = "%s/groundstate_%s.h5" % (groundstateDir, GetGridPrefix(**args))
	filename = GetGroundstateFilename(**args)
	if not os.path.exists(filename):
		raise Exception("Ground state file not found, please run FindGroundstate() (%s)" % filename)
	pyprop.serialization.LoadWavefunctionHDF(filename, "/wavefunction", prop.psi)
	prop.psi.Normalize()

	initPsi = prop.psi.Copy()
	xGrid = prop.psi.GetRepresentation().GetLocalGrid(0)

	normList = []
	corrList = []
	timeList = []

	#handle restarting
	if args.get("restart", False):
		pyprop.PrintOut("Restarting propagation...")
		restartFile = args["restartFile"]

		#Get original propagation end time
		T = prop.Config.Propagation.duration
		
		#load checkpoint wavefunction
		pyprop.serialization.LoadWavefunctionHDF(restartFile, "/wavefunction", prop.psi)

		#get restart time
		restartTime = numpy.zeros(1, dtype=numpy.double)
		if pyprop.ProcId == 0:
			h5file = tables.openFile(restartFile)
			try:
				restartTime[0] = h5file.root.wavefunction._v_attrs.Time
			finally:
				h5file.close()

		#distribute restart time to all procs
		pypar.broadcast(restartTime, 0)

		#restart all propagators with correct time
		prop.RestartPropagation(prop.TimeStep, restartTime[0], T)

	for t in prop.Advance(numOutput):
		N = prop.psi.GetNorm()
		C = abs(prop.psi.InnerProduct(initPsi))**2
		pyprop.PrintOut("t = %f, norm = %f, corr = %f" % (t, N, C))

		normList += [N]
		corrList += [C]
		timeList += [t]

		#other propagation tasks
		for task in propTasks:
			task(t, prop)
	
	prop.NormList = normList
	prop.CorrList = corrList
	prop.TimeList = timeList

	return prop


def Benchmark():
	timer = pyprop.Timers()
	
	timer["Setup"].Start()
	conf = SetupConfig()
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	timer["Setup"].Stop()

	timer["TimeStep"].Start()
	prop.AdvanceStep()
	timer["TimeStep"].Stop()

	timer["TimeStep2"].Start()
	prop.AdvanceStep()
	timer["TimeStep2"].Stop()

	print timer
	
