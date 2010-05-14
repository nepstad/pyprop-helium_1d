import os
import pypar
import pyprop
import tables

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
def GetGroundstateFilename(**args):
	prefix = GetGridPrefix(**args)
	filename = "%s/groundstate_%s.h5" % (GROUNDSTATE_PATH, prefix)
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
		restartTime = -1
		if pyprop.ProcId == 0:
			h5file = tables.openFile(restartFile)
			try:
				restartTime = h5file.root.wavefunction._v_attrs.Time
			finally:
				h5file.close()

		#distribute restart time to all procs
		pypar.broadcast(restartTime, 0)

		#restart all propagators with correct time
		prop.RestartPropagation(prop.TimeStep, restartTime, T)

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
	
