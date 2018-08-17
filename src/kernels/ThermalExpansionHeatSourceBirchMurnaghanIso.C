#include "ThermalExpansionHeatSourceBirchMurnaghanIso.h"

template <>
InputParameters
validParams<ThermalExpansionHeatSourceBirchMurnaghanIso>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Thermal expansion heat source kernel"
                             "Birch Murnaghan equation of state"
                             "Luscher2017 finite strain formalism"
                             "isotropic elastic material");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  return params;
}

ThermalExpansionHeatSourceBirchMurnaghanIso::ThermalExpansionHeatSourceBirchMurnaghanIso(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _thermal_expansion(getParam<Real>("thermal_expansion")), // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _deformation_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient
    _strain_rate(getMaterialProperty<RankTwoTensor>(_base_name + "strain_rate")), // strain rate = dot(epsilon)
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
ThermalExpansionHeatSourceBirchMurnaghanIso::computeQpResidual()
{
  Real heat_source, Je;

  // preliminary calculation of deformation-related quantities
  Je = _deformation_gradient[_qp].det(); // Jacobian = relative volume

  // Equation of state contribution (Psi_EOS in Luscher2017, equation 16)
  // heat source = - K alpha T exp(-n alpha (T - T0)) Je Tr(dot(epsilon))
  heat_source = - _Bulk_Modulus_Ref * _thermal_expansion * _u[_qp]
              * std::exp(-_n_Murnaghan * _thermal_expansion * (_u[_qp] - _reference_temperature))
              * Je * _strain_rate[_qp].trace();

  return - heat_source * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceBirchMurnaghanIso::computeQpJacobian()
{
  // approximate (small strain) Jacobian, first order in the strain rate and temperature
  return _Bulk_Modulus_Ref * _thermal_expansion * _strain_rate[_qp].trace() * _phi[_j][_qp] * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceBirchMurnaghanIso::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val = 0.0;

  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k]) {
      // approximate (small strain) off-diagonal Jacobian, first order in the strain rate and temperature
      val = _Bulk_Modulus_Ref * _thermal_expansion * _u[_qp] * _grad_phi[_j][_qp](k) * _test[_i][_qp] / _dt;
    }
  }

  return val;
}
