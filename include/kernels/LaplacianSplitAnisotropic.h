/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef LAPLACIANSPLITANISOTROPIC_H
#define LAPLACIANSPLITANISOTROPIC_H

#include "KernelGrad.h"

// Forward Declarations
class LaplacianSplitAnisotropic;

template <>
InputParameters validParams<LaplacianSplitAnisotropic>();

/**
 * Split with a variable that holds the Laplacian of the phase field.
 * anisotropic Laplacian for anisotropic damage model
 */
class LaplacianSplitAnisotropic : public KernelGrad
{
public:
  LaplacianSplitAnisotropic(const InputParameters & parameters);

protected:
  virtual RealGradient precomputeQpResidual();
  virtual RealGradient precomputeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _var_c;
  const VariableGradient & _grad_c;
};

#endif // LAPLACIANSPLITANISOTROPIC_H
