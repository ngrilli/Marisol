/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PFTHERMALCONDUCTIVITY_H
#define PFTHERMALCONDUCTIVITY_H

#include "Material.h"

// Forward Declarations
class PFThermalConductivity;
class Function;

template <>
InputParameters validParams<PFThermalConductivity>();

/**
 * Simple material with constant properties.
 */
class PFThermalConductivity : public Material
{
public:
  PFThermalConductivity(const InputParameters & parameters);

protected:
  virtual void computeProperties();

  const Real _k_m;

  const VariableValue & _c;

  const Real _k_c;

  MaterialProperty<Real> & _k;
};

#endif // PFTHERMALCONDUCTIVITY_H
