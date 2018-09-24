//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeFiniteStrainElasticStressGasJWLSolid.h"

registerMooseObject("TensorMechanicsApp", ComputeFiniteStrainElasticStressGasJWLSolid);

template <>
InputParameters
validParams<ComputeFiniteStrainElasticStressGasJWLSolid>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Compute stress using elasticity for finite strains"
                             "third-order Birch Murnaghan equation of state (isotropic)"
                             "for the solid phase"
                             "JWL equation of state"
                             "for the gas phase"
                             "intermediates have gas phase equation of state"
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
  params.addRequiredParam<Real>("A_JWL", "A coefficient in JWL EOS (Lee1973JWL)");
  params.addRequiredParam<Real>("B_JWL", "B coefficient in JWL EOS (Lee1973JWL)");
  params.addRequiredParam<Real>("R1_JWL", "R1 coefficient in JWL EOS (Lee1973JWL)");
  params.addRequiredParam<Real>("R2_JWL", "R2 coefficient in JWL EOS (Lee1973JWL)");
  params.addRequiredParam<Real>("omega_JWL", "omega coefficient in JWL EOS (Lee1973JWL)");
  params.addRequiredParam<Real>("Cv_JWL", "heat capacity coefficient in JWL EOS (Lee1973JWL): E = Cv(T-T0)");
  // other parameters
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  return params;
}

ComputeFiniteStrainElasticStressGasJWLSolid::ComputeFiniteStrainElasticStressGasJWLSolid(
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
    _A_JWL(getParam<Real>("A_JWL")),
    _B_JWL(getParam<Real>("B_JWL")),
    _R1_JWL(getParam<Real>("R1_JWL")),
    _R2_JWL(getParam<Real>("R2_JWL")),
    _omega_JWL(getParam<Real>("omega_JWL")),
    _Cv_JWL(getParam<Real>("Cv_JWL")),
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
ComputeFiniteStrainElasticStressGasJWLSolid::initialSetup()
{
}

void
ComputeFiniteStrainElasticStressGasJWLSolid::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
}

void
ComputeFiniteStrainElasticStressGasJWLSolid::computeQpStress()
{
  Real trD, Kb, KT0prime, Je, V0V, omegaV;
  Real Peos_BM, Peos_JWL; // Birch Murnaghan and JWL pressures
  Real temp = _temp[_qp];
  Real mass_fraction_2_and_3; // mass fraction of final gases and intermediates
  Real mass_fraction_1; // mass fraction of the solid
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

  // Cauchy pressure, JWL (Lee, Finger, Collins, 1973)
  // with E = Cv(T-T0) to include temperature dependence
  omegaV = _omega_JWL * V0V;
  Peos_JWL = _A_JWL * (1.0 - omegaV / _R1_JWL) * std::exp(- _R1_JWL * Je);
  Peos_JWL += _B_JWL * (1.0 - omegaV / _R2_JWL) * std::exp(- _R2_JWL * Je);
  Peos_JWL += omegaV * _Cv_JWL * (temp - _reference_temperature);

  // calculate mass fraction of final gases and intermediates
  mass_fraction_1 = _mass_fraction_1[_qp] ;
  mass_fraction_2_and_3 = 1.0 - mass_fraction_1;

  // volumetric Cauchy stress (equation 17 in Luscher2017)
  iden.zero();
  iden.addIa(1.0);
  _stress[_qp] = (mass_fraction_1 * Peos_BM + mass_fraction_2_and_3 * Peos_JWL) * iden;

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
ComputeFiniteStrainElasticStressGasJWLSolid::rotateQpInitialStress()
{
}
