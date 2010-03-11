from numpy import sin, cos, pi, mod, exp, log
from utils import RegisterProjectNamespace, RegisterAll

@RegisterProjectNamespace
@RegisterAll
def LaserFunction(conf, t):
	E0 = conf.amplitude
	omega = conf.omega
	T = conf.pulse_duration
	#return E0/omega * sin(pi*t/T)**2 * cos(omega * t)
	field = 0.0
	if 0 <= t < T:
		field = E0 * sin(pi*t/T)**2 * cos(omega * t)
	return field


@RegisterProjectNamespace
@RegisterAll
def LaserFunctionGaussian(conf, t): 
    #FWHM duration
    tau = conf.fwhm_pulse_duration

    #total pulse duration
    t0 = conf.pulse_duration

    curField = 0.0
    if 0.0 < t <= 2*t0:
        w = conf.omega
        E0 = conf.amplitude
        curField = E0 * exp(-2 * log(2) * (t-t0)**2 / tau**2) * cos(w * (t-t0))

    return curField
