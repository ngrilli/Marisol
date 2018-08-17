#include "InternalEnergyJumpGasSolid.h"

template <>
InputParameters
validParams<InternalEnergyJumpGasSolid>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Heat source kernel"
                             "for the solid-gas mixture"
                             "to take into account the term"
                             "(e_{BM} - e_{MG}) dot(Y_3)"
                             "representing the jump in internal energy"
                             "when solid and intermediates transform into gas"
                             "must be added to the term ThermalExpansionHeatSourceGasSolid"
                             "to be thermodynamically consistent"
                             "Birch Murnaghan equation of state for the solid phase"
                             "Mie Gruneisen equation of state for the gas phase"
                             "Luscher2017 finite strain formalism"
                             "isotropic elastic material"
                             "pressure average weighted with mass fraction");
  params.addRequiredCoupledVar("mass_fraction_2", "Mass fraction of the second specie");
  // solid parameters
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref_solid", "reference bulk modulus of the solid phase");
  params.addRequiredParam<Real>("bulk_modulus_cor_solid", "KT0 prime correction (Menikoff 2001, Yoo 1998)");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  // gas parameters
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("s_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref_gas", "reference bulk modulus = rho_0 c_0^2 (equation 12 in Menon 2014)");
  params.addRequiredParam<Real>("specific_heat_Gruneisen", "Cv in equation 7 in Menon 2014 to fit gas equation of state");
  // other parameters
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  // the following mass_fraction_rates must refer to Y2_omega2
  // for RDX decomposition model with three mass fractions and two steps reaction
  params.addRequiredParam<MaterialPropertyName>(
      "mass_fraction_rate", "Time rate of the mass fraction due to chemistry");
  params.addParam<MaterialPropertyName>(
      "dmass_fraction_rate_dtemperature", "Material property name with derivative of mass_fraction_rate with temperature");
  params.addParam<MaterialPropertyName>(
      "dmass_fraction_rate_dmass_fraction", "Material property name with derivative of mass_fraction_rate with mass fraction");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  return params;
}

InternalEnergyJumpGasSolid::InternalEnergyJumpGasSolid(const InputParameters & parameters)
  : DerivativeMaterialInterface<HeatSource>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mass_fraction_2(coupledValue("mass_fraction_2")),
    _mass_fraction_2_var(coupled("mass_fraction_2")),
    _mass_fraction_2_name(getVar("mass_fraction_2", 0)->name()),
    // solid parameters
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref_Solid(getParam<Real>("bulk_modulus_ref_solid")),
    _Bulk_Modulus_Cor_Solid(getParam<Real>("bulk_modulus_cor_solid")),
    _thermal_expansion(getParam<Real>("thermal_expansion")), // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
    // gas parameters
    _G_Gruneisen(getParam<Real>("G_Gruneisen")),
    _s_UsUp(getParam<Real>("s_UsUp")),
    _bulk_modulus_ref_gas(getParam<Real>("bulk_modulus_ref_gas")),
    _specific_heat_Gruneisen(getParam<Real>("specific_heat_Gruneisen")),
    // other parameters
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _density(getMaterialProperty<Real>("density_name")),
    // _mass_fraction_rate is dot(Y_3) from arrhenius_material_Y2_omega2 (positive prefactor)
    _mass_fraction_rate(getMaterialProperty<Real>(_base_name + "mass_fraction_rate")),
    _dmass_fraction_rate_dtemperature(
        getMaterialPropertyDerivative<Real>(_base_name + "mass_fraction_rate", _var.name())),
    _dmass_fraction_rate_dmass_fraction(
        getMaterialPropertyDerivative<Real>(_base_name + "mass_fraction_rate", _mass_fraction_2_name)),
    _deformation_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
InternalEnergyJumpGasSolid::computeQpResidual()
{
  Real heat_source, Je, V0V, eta, temp;
  Real Kb, KT0prime; // bulk modulus of the solid and correction
  Real eBM, eMG; // internal energies

  // preliminary calculation of deformation-related quantities
  Je = _deformation_gradient[_qp].det(); // Jacobian = relative volume = v / v_0
  V0V = 1.0 / Je; // inverse relative volume
  eta = 1.0 - Je; // eta = 1 - (v / v0) (Menon, 2014)

  Kb = _Bulk_Modulus_Ref_Solid;
  KT0prime = _Bulk_Modulus_Cor_Solid;
  temp = _u[_qp];

  // Birch-Murnaghan equation of state contribution (Psi_EOS in Luscher2017, equation 16)
  // e_{BM} = Psi_{BM} + T s_{BM}
  eBM = (9.0 / 8.0) * Kb - (9.0 / 16.0) * Kb * (KT0prime - 4.0)
      + (9.0 / 8.0) * Kb * (1.0 - (3.0 / 2.0) * (KT0prime - 4.0)) * std::pow(V0V , 4.0/3.0)
      - (9.0 / 4.0) * Kb * (1.0 - (3.0 / 4.0) * (KT0prime - 4.0)) * std::pow(V0V , 2.0/3.0)
      + (9.0 / 16.0) * Kb * (KT0prime - 4.0) * std::pow(V0V , 2.0)
      + (Kb / _n_Murnaghan) * (Je - 1.0)
      * (std::exp(-_n_Murnaghan * _thermal_expansion * (temp - _reference_temperature))
      * (1.0 + _n_Murnaghan * _thermal_expansion * temp) - 1.0);

  // Mie-Gruneisen equation of state contribution (Psi_EOS in Luscher2017, equation 16)
  // e_{MG} = Psi_{MG} + T s_{MG}
  eMG = _G_Gruneisen * _density[_qp] * _specific_heat_Gruneisen * _reference_temperature * std::log(Je)
      + (_bulk_modulus_ref_gas / (2.0 * (_s_UsUp - 1.0))) * (
        (2.0 * _s_UsUp - 2.0 - _G_Gruneisen) * eta / (_s_UsUp * (1.0 - _s_UsUp * eta))
      + ((_G_Gruneisen * (1.0 - 2.0 * _s_UsUp) + 2.0 * (_s_UsUp - 1.0) * (_s_UsUp - 1.0))
        / (_s_UsUp * _s_UsUp * (_s_UsUp - 1.0)))
        * std::log(1.0 - _s_UsUp * eta)
      + (_G_Gruneisen / (_s_UsUp - 1.0)) * std::log(Je)
      );

  // _mass_fraction_rate is dot(Y_3) from arrhenius_material_Y2_omega2 (positive prefactor)
  // for the three species model of RDX
  heat_source = (eBM - eMG) * _mass_fraction_rate[_qp];

  return - heat_source * _test[_i][_qp];
}

Real
InternalEnergyJumpGasSolid::computeQpJacobian()
{
  Real val;
  Real Je, temp;
  Real Kb, KT0prime; // bulk modulus of the solid and correction
  Real eMGminuseBM; // approximate difference between the internal energies

  // preliminary calculation of deformation-related quantities
  Je = _deformation_gradient[_qp].det(); // Jacobian = relative volume = v / v_0

  Kb = _Bulk_Modulus_Ref_Solid;
  KT0prime = _Bulk_Modulus_Cor_Solid;
  temp = _u[_qp];

  // approximate (small strain and small temperature variation) Jacobian,
  // first order in the strain rate and temperature variation
  // Jacobian = -deBM/dT * mass_fraction_rate + (eMG-eBM) * dmass_fraction_rate_dtemperature
  // deBM/dT = - Kb * (Je - 1.0) * n * alpha^2 * T
  val = Kb * (Je - 1.0) * _n_Murnaghan * _thermal_expansion * _thermal_expansion * temp;
  val *= _mass_fraction_rate[_qp];

  // approximate (small strain and small temperature variation)
  // difference between the internal energies
  eMGminuseBM = (_G_Gruneisen * _density[_qp] * _specific_heat_Gruneisen - _thermal_expansion * Kb)
              * (_reference_temperature - temp) * (Je - 1.0)
              + (_bulk_modulus_ref_gas - Kb) * std::pow(Je - 1.0, 2.0) / 2.0;

  val += eMGminuseBM * _dmass_fraction_rate_dtemperature[_qp];

  return val * _phi[_j][_qp] * _test[_i][_qp];
}

Real
InternalEnergyJumpGasSolid::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  Real Je, temp;
  Real Kb, KT0prime; // bulk modulus of the solid and correction
  Real eMGminuseBM; // approximate difference between the internal energies

  // preliminary calculation of deformation-related quantities
  Je = _deformation_gradient[_qp].det(); // Jacobian = relative volume = v / v_0

  Kb = _Bulk_Modulus_Ref_Solid;
  KT0prime = _Bulk_Modulus_Cor_Solid;
  temp = _u[_qp];

  // approximate (small strain and small temperature variation)
  // difference between the internal energies
  eMGminuseBM = (_G_Gruneisen * _density[_qp] * _specific_heat_Gruneisen - _thermal_expansion * Kb)
              * (_reference_temperature - temp) * (Je - 1.0)
              + (_bulk_modulus_ref_gas - Kb) * std::pow(Je - 1.0, 2.0) / 2.0;

  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k]) {
      // approximate (small strain) off-diagonal Jacobian, first order in the strain rate and temperature variation
      // mass_fraction_rate does not depend on displacement
      // only the difference (e_MG - e_BM) is differentiated
      val = (_G_Gruneisen * _density[_qp] * _specific_heat_Gruneisen - _thermal_expansion * Kb)
          * (_reference_temperature - temp)
          + (_bulk_modulus_ref_gas - Kb)
          * (Je - 1.0);

      val *= _mass_fraction_rate[_qp] * _grad_phi[_j][_qp](k) * _test[_i][_qp];
    }
  }

  if (jvar == _mass_fraction_2_var) {
    val = eMGminuseBM * _dmass_fraction_rate_dmass_fraction[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  }

  return val;
}
