#include <core/common.h>
#include <core/wavefunction.h>
#include <core/representation/coupledspherical/coupledrange.h>
#include <core/transform/spherical/shtools.h>

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_errno.h>


using namespace boost::python;
using namespace blitz;

typedef Array<cplx, 2> MatrixType;
typedef Array<cplx, 1> VectorType;


/*
 * Calculate the projection of the wavefunction on a product of single particle
 * 1D wavefunctions 
 *
 *	<x1(1) x2(2) | psi(1,2) >
 *
 * Remarks:
 * - No symmetrization is made on either psi or the radial functions
 * - Integration weights are assumed to have been applied to psi on beforehand
 *
 * A 2D complex array is returned
 * 		rank0: radial function 1 indices
 * 		rank0: radial function 2 indices
 */
Array<cplx, 2> CalculateProjectionOneParticleStates(MatrixType V1, MatrixType V2, Array<cplx, 2> psiData)
{
	int count1 = V1.extent(1);
	int count2 = V2.extent(1);
	int xcount = V1.extent(0);

	blitz::Array<cplx, 2> proj(count1, count2);
	blitz::Array<cplx, 2> tempProj(xcount, count2);

	//project on V2 states
	tempProj = 0;
	cplx curPsi = 0;
	for (int x1=0; x1<xcount; x1++)
	{
		for (int x2=0; x2<xcount; x2++)
		{
			curPsi = psiData(x1, x2);
			for (int i2=0; i2<count2; i2++)
			{
				tempProj(x1, i2) += conj(V2(x2, i2)) * curPsi;
			}
		}
	}

	//project on V1 states
	cplx a = 0;
	for (int x1=0; x1<xcount; x1++)
	{
		for (int i1=0; i1<count1; i1++)
		{
			a = conj(V1(x1, i1));
			for (int i2=0; i2<count2; i2++)
			{
				proj(i1, i2) += a * tempProj(x1, i2);
			}
		}
	}

	return proj;
}

