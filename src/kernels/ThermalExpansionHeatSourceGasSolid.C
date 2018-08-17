#include "ThermalExpansionHeatSourceGasSolid.h"

template <>
InputParameters
validParams<ThermalExpansionHeatSourceGasSolid>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Thermal expansion heat source kernel"
                             "Birch Murnaghan equation of state for the solid phase"
                             "Mie Gruneisen equation of state for the gas phase"
                             "Luscher2017 finite strain formalism"
                             "isotropic elastic material"
                             "pressure average weighted with mass fraction");
  params.addRequiredCoupledVar("mass_fraction_1", "Mass fraction of the first specie");
  params.addRequiredCoupledVar("mass_fraction_2", "Mass fraction of the second specie");
  // solid parameters
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref_solid", "reference bulk modulus of the solid phase");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  // gas parameters
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("specific_heat_Gruneisen", "Cv in equation 7 in Menon 2014 to fit gas equation of state");
  // other parameters
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  return params;
}

ThermalExpansionHeatSourceGasSolid::ThermalExpansionHeatSourceGasSolid(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mass_fraction_1(coupledValue("mass_fraction_1")),
    _mass_fraction_2(coupledValue("mass_fraction_2")),
    // solid parameters
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref_Solid(getParam<Real>("bulk_modulus_ref_solid")),
    _thermal_expansion(getParam<Real>("thermal_expansion")), // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
    // gas parameters
    _G_Gruneisen(getParam<Real>("G_Gruneisen")),
    _specific_heat_Gruneisen(getParam<Real>("specific_heat_Gruneisen")),
    // other parameters
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _density(getMaterialProperty<Real>("density_name")),
    _deformation_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient
    _strain_rate(getMaterialProperty<RankTwoTensor>(_base_name + "strain_rate")), // strain rate = dot(epsilon)
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
ThermalExpansionHeatSourceGasSolid::computeQpResidual()
{
  Real heat_source, Je;
  Real mass_fraction_3; // mass fraction of final gases
  Real mass_fraction_1_and_2;
  Real strain_rate_trace;

  // preliminary calculation of deformation-related quantities
  Je = _deformation_gradient[_qp].det(); // Jacobian = relative volume
  strain_rate_trace = _strain_rate[_qp].trace();

  // mass fraction of final gases
  mass_fraction_1_and_2 = _mass_fraction_1[_qp] + _mass_fraction_2[_qp];
  mass_fraction_3 = 1.0 - mass_fraction_1_and_2;

  // Birch-Murnaghan equation of state contribution (Psi_EOS in Luscher2017, equation 16)
  // weighted with solid + intermediates mass fraction
  // heat source = - K alpha T exp(-n alpha (T - T0)) Je Tr(dot(epsilon))
  heat_source = - _Bulk_Modulus_Ref_Solid * _thermal_expansion * _u[_qp]
              * std::exp(-_n_Murnaghan * _thermal_expansion * (_u[_qp] - _reference_temperature))
              * Je * strain_rate_trace;
  heat_source *= mass_fraction_1_and_2; // weight with solid + intermediates mass fraction

  // Mie-Gruneisen equation of state contribution (Psi_EOS in Luscher2017, equation 16)
  // weighted with the gas mass fraction
  // heat source = - Gamma rho C_v T Tr(dot(epsilon))
  heat_source -= mass_fraction_3 * _G_Gruneisen * _density[_qp] * _specific_heat_Gruneisen
               * _u[_qp] * strain_rate_trace;

  return - heat_source * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceGasSolid::computeQpJacobian()
{
  Real val;
  Real strain_rate_trace;
  Real mass_fraction_3; // mass fraction of final gases
  Real mass_fraction_1_and_2;

  // preliminary calculation of deformation-related quantities
  strain_rate_trace = _strain_rate[_qp].trace();

  // mass fraction of final gases
  mass_fraction_1_and_2 = _mass_fraction_1[_qp] + _mass_fraction_2[_qp];
  mass_fraction_3 = 1.0 - mass_fraction_1_and_2;

  // approximate (small strain) Jacobian, first order in the strain rate and temperature
  // Birch-Murnaghan contribution
  val = _Bulk_Modulus_Ref_Solid * _thermal_expansion * strain_rate_trace;
  val *= mass_fraction_1_and_2;

  // Mie-Gruneisen contribution
  val += mass_fraction_3 * _G_Gruneisen * _density[_qp] * _specific_heat_Gruneisen * strain_rate_trace;

  return val * _phi[_j][_qp] * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceGasSolid::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  Real strain_rate_trace;
  Real mass_fraction_3; // mass fraction of final gases
  Real mass_fraction_1_and_2;

  // preliminary calculation of deformation-related quantities
  strain_rate_trace = _strain_rate[_qp].trace();

  // mass fraction of final gases
  mass_fraction_1_and_2 = _mass_fraction_1[_qp] + _mass_fraction_2[_qp];
  mass_fraction_3 = 1.0 - mass_fraction_1_and_2;

  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k]) {
      // approximate (small strain) off-diagonal Jacobian, first order in the strain rate and temperature
      // Birch-Murnaghan contribution
      val = _Bulk_Modulus_Ref_Solid * _thermal_expansion * _u[_qp]
          * mass_fraction_1_and_2;

      // Mie-Gruneisen contribution
      val += mass_fraction_3 * _G_Gruneisen * _density[_qp] * _specific_heat_Gruneisen
          * _u[_qp];

      val *= _grad_phi[_j][_qp](k) * _test[_i][_qp] / _dt;
    }
  }

  return val;
}
