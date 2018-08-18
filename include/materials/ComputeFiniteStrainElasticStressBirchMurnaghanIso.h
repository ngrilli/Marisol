//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEFINITESTRAINELASTICSTRESSBIRCHMURNAGHANISO_H
#define COMPUTEFINITESTRAINELASTICSTRESSBIRCHMURNAGHANISO_H

#include "ComputeStressBase.h"
#include "GuaranteeConsumer.h"

class ComputeFiniteStrainElasticStressBirchMurnaghanIso;

template <>
InputParameters validParams<ComputeFiniteStrainElasticStressBirchMurnaghanIso>();

/**
 * ComputeFiniteStrainElasticStressBirchMurnaghanIso computes the stress following elasticity
 * theory for finite strains
 * third-order Birch Murnaghan equation of state
 * finite strain formalism as in Luscher2017 (isotropic)
 * no free energy calculation for efficiency
 */
class ComputeFiniteStrainElasticStressBirchMurnaghanIso : public ComputeStressBase, public GuaranteeConsumer
{
public:
  ComputeFiniteStrainElasticStressBirchMurnaghanIso(const InputParameters & parameters);

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

  // exponent in Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus
  const Real _Bulk_Modulus_Ref;

  // KT0 prime correction to the bulk modulus (Menikoff 2001)
  const Real _Bulk_Modulus_Cor;

  // Von Neumann coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

  // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
  const Real _thermal_expansion;

  // reference temperature, as in Luscher2017
  const Real _reference_temperature;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
};

#endif // COMPUTEFINITESTRAINELASTICSTRESSBIRCHMURNAGHANISO_H
