import sys
import os

sys.path.append("./pyprop")
import pyprop
pyprop.ProjectNamespace = globals()

from libpotential import *

from numpy import sin, cos, pi


def SetupConfig(**args):
	configFile = 'config.ini'
	if 'config' in args:
		configFile = args['config']
	conf = pyprop.Load(configFile)

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


def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	prop.psi.Normalize()

	return prop


def GetGridPrefix(**args):
	conf = SetupConfig(**args)
	xmax = conf.Representation.rank0[1]
	count = conf.Representation.rank0[2]
	return "xmax%s_n%s" % (xmax, count)


def FindGroundstate(**args):
	"""
	Use imaginary time propagation to find ground state
	"""
	args["imtime"] = True

	filename = "groundstate/groundstate_%s.h5" % GetGridPrefix(**args)
	dir = os.path.dirname(filename)
	if not os.path.exists(dir) and pyprop.ProcId == 0:
		os.makedirs(dir)

	#Setup problem
	prop = SetupProblem(**args)

	#Propagate in imaginary time to obtain ground state
	for t in prop.Advance(10):
		N = prop.psi.GetNorm()
		E = prop.GetEnergy()
		print "t = %03f, N = %s, E = %s" % (t, N, E)

	#Save groundstate
	prop.psi.Normalize()
	prop.SaveWavefunctionHDF(filename, "/wavefunction")

	return prop


def Propagate(**args):
	args['imtime'] = False
	args['additional_potentials'] = ["LaserPotential"]
	args['config'] = "propagation.ini"
	numOutput = args.get("numOutput", 10)
	
	#Setup problem
	prop = SetupProblem(**args)

	#Get initial state
	filename = "groundstate/groundstate_%s.h5" % GetGridPrefix(**args)
	if not os.path.exists(filename):
		raise Exception("Ground state file not found, please run FindGroundstate() (%s)" % filename)
	pyprop.serialization.LoadWavefunctionHDF(filename, "/wavefunction", prop.psi)
	prop.psi.Normalize()

	initPsi = prop.psi.Copy()

	for t in prop.Advance(numOutput):
		N = prop.psi.GetNorm()
		C = abs(prop.psi.InnerProduct(initPsi))**2
		print "t = %f, norm = %f, corr = %f" % (t, N, C)

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
	

def LaserFunction(conf, t):
	E0 = conf.e0
	omega = conf.omega
	T = conf.pulse_duration
	#return E0/omega * sin(pi*t/T)**2 * cos(omega * t)
	return E0 * sin(pi*t/T)**2 * cos(omega * t)


def LaserFunctionGaussian(conf, t):
	pass
