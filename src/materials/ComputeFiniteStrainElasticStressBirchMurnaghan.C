//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeFiniteStrainElasticStressBirchMurnaghan.h"

registerMooseObject("TensorMechanicsApp", ComputeFiniteStrainElasticStressBirchMurnaghan);

template <>
InputParameters
validParams<ComputeFiniteStrainElasticStressBirchMurnaghan>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Compute stress using elasticity for finite strains"
                             "third-order Birch Murnaghan equation of state"
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

ComputeFiniteStrainElasticStressBirchMurnaghan::ComputeFiniteStrainElasticStressBirchMurnaghan(
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
    _pk2(declareProperty<RankTwoTensor>("pk2")), // 2nd Piola Kirchoff Stress
    _lag_e(declareProperty<RankTwoTensor>("lage")), // Lagrangian strain
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient"))
{
}

void
ComputeFiniteStrainElasticStressBirchMurnaghan::initialSetup()
{
}

void
ComputeFiniteStrainElasticStressBirchMurnaghan::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
}

void
ComputeFiniteStrainElasticStressBirchMurnaghan::computeQpStress()
{
  Real trD, Kb, KT0prime, Je, Peos, V0V, delta;
  Real temp = _temp[_qp];
  RankTwoTensor iden, ce, invce, thermal_eigenstrain;

  // reference bulk modulus is an input
  Kb = _Bulk_Modulus_Ref;
  KT0prime = _Bulk_Modulus_Cor;

  // tensor calculation
  iden.zero();
  iden.addIa(1.0);

  ce = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp]; // Cauchy-Green strain tensor

  _lag_e[_qp] = ce - iden; // Lagrangian strain tensor
  _lag_e[_qp] = _lag_e[_qp] * 0.5;

  invce = ce.inverse();

  Je = _deformation_gradient[_qp].det(); // Jacobian = relative volume

  // Cauchy pressure, third-order Birch Murnaghan (Menikoff2001, Yoo1998)
  V0V = 1.0 / Je; // relative volume
  Peos = 1.5 * Kb * (std::pow(V0V , 7.0/3.0) - std::pow(V0V , 5.0/3.0));
  Peos = Peos * (1.0 + 0.75 * (KT0prime - 4.0) * (std::pow(V0V , 2.0/3.0) - 1.0));
  Peos = - Peos; // negative stress in compression
  Peos += (Kb / _n_Murnaghan) * (std::exp(-_n_Murnaghan * _thermal_expansion * (temp - _reference_temperature)) - 1.0);

  // volumetric stress (equation 17 in Luscher2017)
  _pk2[_qp] = Je * Peos * invce;

  // thermal eigenstrain (equation (18) in Luscher2017)
  // Lagrangian strain E_thermal = 1/2 (F_thermal^T F_thermal - I)
  // finite strain formula (Lubarda2002): F_thermal = exp((alpha/3)*(T-T_ref)) I
  thermal_eigenstrain = (1.0 / 2.0)
                      * (std::exp((2.0/3.0) * _thermal_expansion * (temp - _reference_temperature)) - 1.0)
                      * iden;

  // deviatoric/isochoric stress (equation (18) in Luscher2017): C : (Ee - alpha)
  _pk2[_qp] += _elasticity_tensor[_qp] * (_lag_e[_qp] - thermal_eigenstrain);

  // Pcor = correcting pressure = linearized form of the EOS
  // equation (18) in Luscher2017
  delta = 1.5 * (std::pow(Je , 2.0/3.0) - 1.0);
  _pk2[_qp] -= Kb * std::pow(Je , 2.0/3.0)
           * (delta * iden - 3.0 * thermal_eigenstrain)
           * invce;

  // Calculate bulk viscosity damping
  // C0 * dot(J) / J * |dot(J) / J| + C1 * dot(J) / J
  // C0 should be chosen of the order of rho * Le^2, rho = density, Le = element size
  // C1 should be chosen of the order of rho * Le * cs, cs = sound speed
  // Maheo et al. Mechanics Research Communications 38 (2011) 81 88
  trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  trD /= _deformation_gradient_old[_qp].det();

  _pk2[_qp].addIa( _C0 * trD * std::abs(trD) );
  _pk2[_qp].addIa( _C1 * trD );

  // Cauchy stress
  _stress[_qp] = _deformation_gradient[_qp] * _pk2[_qp] * _deformation_gradient[_qp].transpose() / Je;

  // Compute dstress_dstrain
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp]; // This is NOT the exact jacobian
}

void
ComputeFiniteStrainElasticStressBirchMurnaghan::rotateQpInitialStress()
{
}
