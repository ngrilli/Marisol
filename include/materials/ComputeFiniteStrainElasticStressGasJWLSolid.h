//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEFINITESTRAINELASTICSTRESSGASJWLSOLID_H
#define COMPUTEFINITESTRAINELASTICSTRESSGASJWLSOLID_H

#include "ComputeStressBase.h"
#include "GuaranteeConsumer.h"

class ComputeFiniteStrainElasticStressGasJWLSolid;

template <>
InputParameters validParams<ComputeFiniteStrainElasticStressGasJWLSolid>();

/**
 * ComputeFiniteStrainElasticStressGasJWLSolid computes the stress following elasticity
 * theory for finite strains
 * third-order Birch Murnaghan equation of state (isotropic)
 * for the solid phase
 * JWL equation of state
 * (Lee, Finger, Collins, JWL equation of state coefficients for high explosives, 1973)
 * for the gas phase
 * intermediates have gas phase equation of state
 * finite strain formalism as in Luscher2017
 * no free energy calculation for efficiency
 */
class ComputeFiniteStrainElasticStressGasJWLSolid : public ComputeStressBase, public GuaranteeConsumer
{
public:
  ComputeFiniteStrainElasticStressGasJWLSolid(const InputParameters & parameters);

  void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

  /**
   * InitialStress Deprecation: remove this method
   *
   * Rotates initial_stress via rotation_increment.
   * In large-strain scenarios this must be used before addQpInitialStress
   */
  virtual void rotateQpInitialStress();

  // temperature
  const VariableValue & _temp;

  // Mass fraction of the first specie
  const VariableValue & _mass_fraction_1;

  // Mass fraction of the second specie
  const VariableValue & _mass_fraction_2;

  // exponent in Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus of the solid
  const Real _Bulk_Modulus_Ref_Solid;

  // KT0 prime correction to the bulk modulus of the solid (Menikoff 2001)
  const Real _Bulk_Modulus_Cor_Solid;

  // volumetric thermal exapnsion coeffcient of the solid, as in Austin Barton 2015
  const Real _thermal_expansion;

  // A coefficient in JWL EOS (Lee1973JWL)
  const Real _A_JWL;

  // B coefficient in JWL EOS (Lee1973JWL)
  const Real _B_JWL;

  // R1 coefficient in JWL EOS (Lee1973JWL)
  const Real _R1_JWL;

  // R2 coefficient in JWL EOS (Lee1973JWL)
  const Real _R2_JWL;

  // omega coefficient in JWL EOS (Lee1973JWL)
  const Real _omega_JWL;

  // heat capacity coefficient in JWL EOS (Lee1973JWL): E = Cv(T-T0)
  const Real _Cv_JWL;

  // reference temperature, as in Luscher2017
  const Real _reference_temperature;

  // Von Neumann coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
};

#endif // COMPUTEFINITESTRAINELASTICSTRESSGASJWLSOLID_H
