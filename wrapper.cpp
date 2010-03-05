
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <analysis.cpp>
#include <potential.cpp>
#include <potential1d.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpotential)
{
    class_< DynamicPotentialEvaluator<SoftColoumbPotential<2>,2>, boost::noncopyable >("SoftColoumbPotential_2", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<SoftColoumbPotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<SoftColoumbPotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<SoftColoumbPotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<SoftColoumbPotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<SoftColoumbPotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<SoftColoumbPotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<SoftColoumbPotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<LaserDipolePotential<2>,2>, boost::noncopyable >("LaserDipolePotential_2", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<LaserDipolePotential<2>,2>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<LaserDipolePotential<2>,2>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<LaserDipolePotential<2>,2>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<LaserDipolePotential<2>,2>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<LaserDipolePotential<2>,2>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<LaserDipolePotential<2>,2>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<LaserDipolePotential<2>,2>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<LaserDipolePotential<1>,1>, boost::noncopyable >("LaserDipolePotential_1", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<LaserDipolePotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<LaserDipolePotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<LaserDipolePotential<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<LaserDipolePotential<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<LaserDipolePotential<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<LaserDipolePotential<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<LaserDipolePotential<1>,1>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<OneElectronSoftColoumbPotential<1>,1>, boost::noncopyable >("OneElectronSoftColoumbPotential_1", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<OneElectronSoftColoumbPotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<OneElectronSoftColoumbPotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<OneElectronSoftColoumbPotential<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<OneElectronSoftColoumbPotential<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<OneElectronSoftColoumbPotential<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<OneElectronSoftColoumbPotential<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<OneElectronSoftColoumbPotential<1>,1>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<ComplexAbsorbingPotential<1>,1>, boost::noncopyable >("ComplexAbsorbingPotential_1", init<  >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<ComplexAbsorbingPotential<1>,1>::CalculateExpectationValue)
    ;

    def("CalculateProjectionOneParticleStates", CalculateProjectionOneParticleStates);
}

