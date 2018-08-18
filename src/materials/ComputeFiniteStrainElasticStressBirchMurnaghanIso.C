//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeFiniteStrainElasticStressBirchMurnaghanIso.h"

registerMooseObject("TensorMechanicsApp", ComputeFiniteStrainElasticStressBirchMurnaghanIso);

template <>
InputParameters
validParams<ComputeFiniteStrainElasticStressBirchMurnaghanIso>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Compute stress using elasticity for finite strains"
                             "third-order Birch Murnaghan equation of state (isotropic)"
                             "finite strain formalism as in Luscher2017"
                             "no free energy calculation for efficiency");
  params.addRequiredCoupledVar("temp","Temperature");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("bulk_modulus_cor", "KT0 prime correction (Menikoff 2001, Yoo 1998)");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  return params;
}

ComputeFiniteStrainElasticStressBirchMurnaghanIso::ComputeFiniteStrainElasticStressBirchMurnaghanIso(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
    GuaranteeConsumer(this),
    _temp(coupledValue("temp")),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _Bulk_Modulus_Cor(getParam<Real>("bulk_modulus_cor")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _thermal_expansion(getParam<Real>("thermal_expansion")), // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient"))
{
}

void
ComputeFiniteStrainElasticStressBirchMurnaghanIso::initialSetup()
{
}

void
ComputeFiniteStrainElasticStressBirchMurnaghanIso::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
}

void
ComputeFiniteStrainElasticStressBirchMurnaghanIso::computeQpStress()
{
  Real trD, Kb, KT0prime, Je, Peos, V0V;
  Real temp = _temp[_qp];
  RankTwoTensor iden; // identity matrix

  // reference bulk modulus is an input
  Kb = _Bulk_Modulus_Ref;
  KT0prime = _Bulk_Modulus_Cor;

  Je = _deformation_gradient[_qp].det(); // Jacobian = relative volume

  // Cauchy pressure, third-order Birch Murnaghan (Menikoff2001, Yoo1998)
  V0V = 1.0 / Je; // relative volume
  Peos = 1.5 * Kb * (std::pow(V0V , 7.0/3.0) - std::pow(V0V , 5.0/3.0));
  Peos = Peos * (1.0 + 0.75 * (KT0prime - 4.0) * (std::pow(V0V , 2.0/3.0) - 1.0));
  Peos = - Peos; // negative stress in compression
  Peos += (Kb / _n_Murnaghan) * (std::exp(-_n_Murnaghan * _thermal_expansion * (temp - _reference_temperature)) - 1.0);

  iden.zero();
  iden.addIa(1.0);

  // volumetric Cauchy stress (equation 17 in Luscher2017)
  _stress[_qp] = Peos * iden;

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
ComputeFiniteStrainElasticStressBirchMurnaghanIso::rotateQpInitialStress()
{
}
