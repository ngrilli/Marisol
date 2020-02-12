//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ACGRGRPLASTICDRIVINGFORCE_H
#define ACGRGRPLASTICDRIVINGFORCE_H

#include "ACBulk.h"

// Forward Declarations
class ACGrGrPlasticDrivingForce;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
template <typename>
class RankFourTensorTempl;
typedef RankFourTensorTempl<Real> RankFourTensor;

//class RankTwoTensor;
//class RankFourTensor;

template <>
InputParameters validParams<ACGrGrPlasticDrivingForce>();

/**
 * Calculates the porton of the Allen-Cahn equation that results from the deformation energy.
 * Must access the elastic_strain stored as a material property
 * Requires the name of the elastic tensor derivative and the plastic energy as an input.
 */
class ACGrGrPlasticDrivingForce : public ACBulk<Real>
{
public:
  ACGrGrPlasticDrivingForce(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);

private:
  const MaterialProperty<Real> & _D_plastic_energy;
  const MaterialProperty<RankFourTensor> & _D_elastic_tensor;
  const MaterialProperty<RankTwoTensor> & _elastic_strain;
};

#endif // ACGRGRPLASTICDRIVINGFORCE_H
