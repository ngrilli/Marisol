/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputeThermalExpansionEigenstrainAnisotropic.h"
#include "RotationTensor.h"

template<>
InputParameters validParams<ComputeThermalExpansionEigenstrainAnisotropic>()
{
  InputParameters params = validParams<ComputeThermalExpansionEigenstrainAnisotropicBase>();
  params.addClassDescription("Computes eigenstrain due to thermal expansion with a constant coefficient");
  params.addParam<std::vector<Real>>("thermal_expansion_coeff", "Thermal expansion coefficient");

  return params;
}

ComputeThermalExpansionEigenstrainAnisotropic::ComputeThermalExpansionEigenstrainAnisotropic(const InputParameters & parameters) :
    ComputeThermalExpansionEigenstrainAnisotropicBase(parameters),
    _thermal_expansion_coeff(getParam<std::vector<Real>>("thermal_expansion_coeff"))
{
}


void
ComputeThermalExpansionEigenstrainAnisotropic::computeThermalStrain(std::vector<Real> & thermal_strain, std::vector<Real> & instantaneous_cte)
{
  thermal_strain[0] = _thermal_expansion_coeff[0] * (_temperature[_qp] - _stress_free_temperature);
  thermal_strain[1] = _thermal_expansion_coeff[1] * (_temperature[_qp] - _stress_free_temperature);
  thermal_strain[2] = _thermal_expansion_coeff[2] * (_temperature[_qp] - _stress_free_temperature);
  
  instantaneous_cte[0] = _thermal_expansion_coeff[0];
  instantaneous_cte[1] = _thermal_expansion_coeff[1];
  instantaneous_cte[2] = _thermal_expansion_coeff[2];
}



