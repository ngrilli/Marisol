//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeFiniteStrainElasticStressGasSolid.h"

registerMooseObject("TensorMechanicsApp", ComputeFiniteStrainElasticStressGasSolid);

template <>
InputParameters
validParams<ComputeFiniteStrainElasticStressGasSolid>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Compute stress using elasticity for finite strains"
                             "third-order Birch Murnaghan equation of state (isotropic)"
                             "for the solid phase"
                             "Mie Gruneisen equation of state"
                             "for the gas phase"
                             "finite strain formalism as in Luscher2017"
                             "no free energy calculation for efficiency");
  params.addRequiredCoupledVar("temp","Temperature");
  params.addRequiredCoupledVar("mass_fraction_1", "Mass fraction of the first specie");
  params.addRequiredCoupledVar("mass_fraction_2", "Mass fraction of the second specie");
  // solid parameters
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref_solid", "reference bulk modulus");
  params.addRequiredParam<Real>("bulk_modulus_cor_solid", "KT0 prime correction (Menikoff 2001, Yoo 1998)");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  // gas parameters
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("s_UsUp", "Us-Up slope in Mie-Gruneisen EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref_gas", "reference bulk modulus = rho_0 c_0^2 (equation 12 in Menon 2014)");
  params.addRequiredParam<Real>("specific_heat_Gruneisen", "Cv in equation 7 in Menon 2014 to fit gas equation of state");
  // other parameters
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  return params;
}

ComputeFiniteStrainElasticStressGasSolid::ComputeFiniteStrainElasticStressGasSolid(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
    GuaranteeConsumer(this),
    _temp(coupledValue("temp")),
    _mass_fraction_1(coupledValue("mass_fraction_1")),
    _mass_fraction_2(coupledValue("mass_fraction_2")),
    // solid parameters
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref_Solid(getParam<Real>("bulk_modulus_ref_solid")),
    _Bulk_Modulus_Cor_Solid(getParam<Real>("bulk_modulus_cor_solid")),
    _thermal_expansion(getParam<Real>("thermal_expansion")), // volumetric thermal expansion coeffcient, as in Austin Barton 2015
    // gas parameters
    _G_Gruneisen(getParam<Real>("G_Gruneisen")),
    _s_UsUp(getParam<Real>("s_UsUp")),
    _bulk_modulus_ref_gas(getParam<Real>("bulk_modulus_ref_gas")),
    _specific_heat_Gruneisen(getParam<Real>("specific_heat_Gruneisen")),
    // other parameters
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _density(getMaterialProperty<Real>("density_name")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient"))
{
}

void
ComputeFiniteStrainElasticStressGasSolid::initialSetup()
{
}

void
ComputeFiniteStrainElasticStressGasSolid::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
}

void
ComputeFiniteStrainElasticStressGasSolid::computeQpStress()
{
  Real trD, Kb, KT0prime, Je, V0V, eta;
  Real Peos_BM, Peos_MG; // Birch Murnaghan and Mie Gruneisen pressures
  Real temp = _temp[_qp];
  Real mass_fraction_3; // mass fraction of final gases
  Real mass_fraction_1_and_2;
  RankTwoTensor iden; // identity matrix

  // reference bulk modulus is an input
  Kb = _Bulk_Modulus_Ref_Solid;
  KT0prime = _Bulk_Modulus_Cor_Solid;

  Je = _deformation_gradient[_qp].det(); // Jacobian = relative volume

  // Cauchy pressure, third-order Birch Murnaghan (Menikoff2001, Yoo1998)
  V0V = 1.0 / Je; // relative volume
  Peos_BM = 1.5 * Kb * (std::pow(V0V , 7.0/3.0) - std::pow(V0V , 5.0/3.0));
  Peos_BM = Peos_BM * (1.0 + 0.75 * (KT0prime - 4.0) * (std::pow(V0V , 2.0/3.0) - 1.0));
  Peos_BM = - Peos_BM; // negative stress in compression
  Peos_BM += (Kb / _n_Murnaghan) * (std::exp(-_n_Murnaghan * _thermal_expansion * (temp - _reference_temperature)) - 1.0);

  // Cauchy pressure, Mie Gruneisen (Akiki, Menon, 2017)
  eta = 1.0 - Je; // eta = 1 - (v / v0) (Menon, 2014)
  Peos_MG = _G_Gruneisen * _density[_qp] * _specific_heat_Gruneisen * (temp - _reference_temperature) * V0V;
  Peos_MG += _bulk_modulus_ref_gas * eta * (1.0 - (_G_Gruneisen / 2.0) * (V0V - 1.0)) / std::pow((1.0 - _s_UsUp * eta) , 2.0);
  Peos_MG = - Peos_MG; // negative stress in compression

  // mass fraction of final gases
  mass_fraction_1_and_2 = _mass_fraction_1[_qp] + _mass_fraction_2[_qp];
  mass_fraction_3 = 1.0 - mass_fraction_1_and_2;

  // volumetric Cauchy stress (equation 17 in Luscher2017)
  iden.zero();
  iden.addIa(1.0);
  _stress[_qp] = (mass_fraction_1_and_2 * Peos_BM + mass_fraction_3 * Peos_MG) * iden;

  // Calculate bulk viscosity damping
  // C0 * dot(J) / J * |dot(J) / J| + C1 * dot(J) / J
  // C0 should be chosen of the order of rho * Le^2, rho = density, Le = element size
  // C1 should be chosen of the order of rho * Le * cs, cs = sound speed
  // Maheo et al. Mechanics Research Communications 38 (2011) 81 88
  trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  trD /= _deformation_gradient_old[_qp].det();

  _stress[_qp].addIa( _C0 * trD * std::abs(trD) );
  _stress[_qp].addIa( _C1 * trD );

  // Compute dstress_dstrain
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp]; // This is NOT the exact jacobian
}

void
ComputeFiniteStrainElasticStressGasSolid::rotateQpInitialStress()
{
}
