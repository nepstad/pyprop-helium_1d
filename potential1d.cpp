#include <core/wavefunction.h>
#include <core/potential/dynamicpotentialevaluator.h>

template<int Rank>
class OneElectronSoftColoumbPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	//Potential parameters
	double nuclearSoft;
	double nuclearCharge;
	double electronSoft;
	double electronCharge;

	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("nuclear_soft", nuclearSoft);
		config.Get("nuclear_charge", nuclearCharge);
	}

	inline double GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = pos(0);

		double V = nuclearCharge / sqrt(x*x + nuclearSoft);

		return V;
	}
};


template<int Rank>
class ComplexAbsorbingPotential : public PotentialBase<Rank>
{
public:
	//Required by DynamicPotentialEvaluator
	cplx TimeStep;
	double CurTime;

	double scalingReal;
	double scalingImag;
	double factorReal;
	double factorImag;
	double absorberStart;
	double absorberLength;

	/*
	 * Called once with the corresponding config section
	 * from the configuration file. Do all one time set up routines
	 * here.
	 */
	void ApplyConfigSection(const ConfigSection &config)
	{
		config.Get("absorber_start", absorberStart);
		config.Get("absorber_length", absorberLength);
		config.Get("scaling_real", scalingReal);
		config.Get("scaling_imag", scalingImag);
		config.Get("factor_real", factorReal);
		config.Get("factor_imag", factorImag);
	}

	/*
	 * Called for every grid point 
	 */
	inline cplx GetPotentialValue(const blitz::TinyVector<double, Rank> &pos)
	{
		double x = abs(pos(0));
		cplx V = 0;
		if (x > absorberStart)
		{
			double curLength = (x - absorberStart) / absorberLength;
			double Vr = factorReal * std::pow(curLength, scalingReal);
			double Vi = factorImag * std::pow(curLength, scalingImag);
			V += cplx(Vr , Vi);
		}
		return V;
	}
};

