/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "LaplacianSplitAnisotropic.h"

template <>
InputParameters
validParams<LaplacianSplitAnisotropic>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Split with a variable that holds the Laplacian of a phase field variable."
      "anisotropic Laplacian for anisotropic damage model (Nguyen, 2017)");
  params.addRequiredCoupledVar("c", "Field variable to take the Laplacian of");
  params.addRequiredParam<Real>("beta_penalty","penalty for damage on planes not normal to z (Nguyen, 2017)");
  params.addRequiredParam<std::vector<Real>>("cleavage_plane_normal","Normal to the favoured cleavage plane");
  return params;
}

LaplacianSplitAnisotropic::LaplacianSplitAnisotropic(const InputParameters & parameters)
  : Kernel(parameters),
    _var_c(coupled("c")),
    _grad_c(coupledGradient("c")),
    _beta_penalty(getParam<Real>("beta_penalty")),
    _cleavage_plane_normal(getParam<std::vector<Real>>("cleavage_plane_normal"))
{
}

Real
LaplacianSplitAnisotropic::computeQpResidual()
{
  Real Mcoupling; // scalar product: M_i grad_i c

  Mcoupling = _grad_c[_qp](0) * _cleavage_plane_normal[0]
            + _grad_c[_qp](1) * _cleavage_plane_normal[1]
            + _grad_c[_qp](2) * _cleavage_plane_normal[2];

  return (1.0 + _beta_penalty) * _grad_c[_qp] * _grad_test[_i][_qp]
         - _beta_penalty * Mcoupling * (_cleavage_plane_normal[0] * _grad_test[_i][_qp](0)
                                      + _cleavage_plane_normal[1] * _grad_test[_i][_qp](1)
                                      + _cleavage_plane_normal[2] * _grad_test[_i][_qp](2));
}

Real
LaplacianSplitAnisotropic::computeQpJacobian()
{
  return 0.0;
}

Real
LaplacianSplitAnisotropic::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _var_c)
    return _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  return 0.0;
}
