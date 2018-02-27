/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef LAPLACIANSPLITANISOTROPIC_H
#define LAPLACIANSPLITANISOTROPIC_H

#include "Kernel.h"

// Forward Declarations
class LaplacianSplitAnisotropic;

template <>
InputParameters validParams<LaplacianSplitAnisotropic>();

/**
 * Split with a variable that holds the Laplacian of the phase field.
 * anisotropic Laplacian for anisotropic damage model
 */
class LaplacianSplitAnisotropic : public Kernel
{
public:
  LaplacianSplitAnisotropic(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const unsigned int _var_c;
  const VariableGradient & _grad_c;

  // penalty for damage on planes not normal to the favoured cleavage plane normal (Nguyen, 2017)
  const Real _beta_penalty;

  // normal to the favoured cleavage plane: M in (Nguyen, 2017)
  const std::vector<Real> _cleavage_plane_normal;
};

#endif // LAPLACIANSPLITANISOTROPIC_H
