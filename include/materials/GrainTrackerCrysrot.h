/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef GRAINTRACKERCRYSROT_H
#define GRAINTRACKERCRYSROT_H

#include "GrainDataTracker.h"
#include "EulerAngleProvider.h"
#include "RankTwoTensor.h"


class GrainTrackerCrysrot;
class EulerAngleProvider;

template<>
InputParameters validParams<GrainTrackerCrysrot>();

/**
 * Manage a list of elasticity tensors for the grains
 */
class GrainTrackerCrysrot : public GrainDataTracker<EulerAngles>
{
public:
  GrainTrackerCrysrot(const InputParameters & parameters);

protected:
  EulerAngles newGrain(unsigned int new_grain_id);

  /// generate random rotations when the Euler Angle provider runs out of data (otherwise error out)
  const bool _random_rotations;

  /// unrotated elasticity tensor

  /// object providing the Euler angles
  const EulerAngleProvider & _euler;
};

#endif // GRAINTRACKERCRYSROT_H
