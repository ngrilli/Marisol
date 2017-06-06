/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PFThermalConductivity.h"

// libmesh includes
#include "libmesh/quadrature.h"

template <>
InputParameters validParams<PFThermalConductivity>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("material_thermal_conductivity",
                                "Thermal conductivity of the undamaged material");
  params.addRequiredCoupledVar("c", "Coupled phase field");

  params.addParam<Real>("crack_thermal_conductivity", 0.0,
                        "Thermal conductivity of fully damaged material");
  params.addClassDescription("This material calculates effective thermal"
                             "conductivity of partially damaged material"
                             "k_eff = (1-c^2)*k_m + c^2*k_c");
  return params;
}

PFThermalConductivity::PFThermalConductivity(
    const InputParameters & parameters)
  : Material(parameters),

    _k_m(getParam<Real>("material_thermal_conductivity")),
    _c(coupledValue("c")),
    _k_c(getParam<Real>("crack_thermal_conductivity")),
    _k(declareProperty<Real>("thermal_conductivity"))
{
}

void PFThermalConductivity::computeProperties()
{
  for (unsigned int qp(0); qp < _qrule->n_points(); ++qp)
  {
    Real qp_c = _c[qp];
      if (_c[qp] < 0) { qp_c = 0;}
      if (_c[qp] > 1) {qp_c = 1;}

    _k[qp] = _k_m*(1-(qp_c*qp_c)) + _k_c*qp_c*qp_c;
  }
}
