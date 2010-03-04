from numpy import sin, cos, pi, mod
from utils import RegisterProjectNamespace, RegisterAll

@RegisterProjectNamespace
@RegisterAll
def LaserFunction(conf, t):
	E0 = conf.e0
	omega = conf.omega
	T = conf.pulse_duration
	#return E0/omega * sin(pi*t/T)**2 * cos(omega * t)
	field = 0.0
	if 0 <= t < T:
		field = E0 * sin(pi*t/T)**2 * cos(omega * t)
	return field


def LaserFunctionGaussian(conf, t):
	pass
