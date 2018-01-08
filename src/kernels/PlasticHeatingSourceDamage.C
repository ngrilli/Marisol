#include "PlasticHeatingSourceDamage.h"

template <>
InputParameters
validParams<PlasticHeatingSourceDamage>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Heat source kernel from plasticity with damage (Miehe 2016)");
  params.addRequiredParam<Real>("plastic_factor","Prefactor of the plastic contribution to damage");
  params.addParam<MaterialPropertyName>("W0p","plastic energy (unbroken)");
  params.addParam<MaterialPropertyName>("W0p_broken","plastic energy (broken)");
  params.addParam<MaterialPropertyName>("dW0p_dstrain", "Derivative of W0p (unbroken) with strain");
  params.addParam<MaterialPropertyName>("dW0p_broken_dstrain", "Derivative of W0p (broken) with strain");
  params.addRequiredCoupledVar("c", "Damage phase field");
  params.addCoupledVar("displacements",
                       "The string of displacements suitable for the problem statement");
  return params;
}

PlasticHeatingSourceDamage::PlasticHeatingSourceDamage(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _plastic_factor(getParam<Real>("plastic_factor")), // Pre-Factor of the plastic energy contribution to damage evolution
    _W0p(getMaterialProperty<Real>("W0p")), // plastic energy (unbroken)
    _W0p_old(getMaterialPropertyOld<Real>("W0p")), // old plastic energy (unbroken)
    _W0p_broken(getMaterialProperty<Real>("W0p_broken")), // plastic energy (broken)
    _W0p_broken_old(getMaterialPropertyOld<Real>("W0p_broken")), // old plastic energy (broken)
    _dW0p_dstrain(isParamValid("dW0p_dstrain") ? &getMaterialProperty<RankTwoTensor>("dW0p_dstrain"): NULL),
    _dW0p_broken_dstrain(isParamValid("dW0p_broken_dstrain") ? &getMaterialProperty<RankTwoTensor>("dW0p_broken_dstrain"): NULL),
    _c(coupledValue("c")), // order parameter for damage
    _c_var(coupled("c")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
PlasticHeatingSourceDamage::computeQpResidual()
{
  Real heat_source = 0.0;

  heat_source = _plastic_factor * (1.0 - _c[_qp]) * (1.0 - _c[_qp])
                * ( _W0p[_qp] - _W0p_old[_qp] ) / _dt * _test[_i][_qp]; // release of stored energy due to plastic deformation
  heat_source += ( _W0p_broken[_qp] - _W0p_broken_old[_qp] ) / _dt * _test[_i][_qp]; // plastic work (broken)

  return - heat_source;
}

Real
PlasticHeatingSourceDamage::computeQpJacobian()
{
  return 0.0 * _phi[_j][_qp] * _test[_i][_qp];
}

Real
PlasticHeatingSourceDamage::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val = 0.0;
  unsigned int c_comp;
  bool disp_flag = false;

  val = 0.0;
  // Contribution of auxiliary variable to off diag Jacobian of temperature kernel
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k])
    {
      c_comp = k;
      disp_flag = true;
    }
  }

  // Contribution of displacements to off diag Jacobian of temperature kernel (unbroken)
  if (disp_flag && _dW0p_dstrain != NULL)
  {
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      val += _plastic_factor * (1.0 - _c[_qp]) * (1.0 - _c[_qp])
             * ( (*_dW0p_dstrain)[_qp](c_comp, i) + (*_dW0p_dstrain)[_qp](i, c_comp) ) * _grad_phi[_j][_qp](i) / _dt;
    }
  }

  // Contribution of displacements to off diag Jacobian of temperature kernel (broken)
  if (disp_flag && _dW0p_broken_dstrain != NULL)
  {
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      val += ( (*_dW0p_broken_dstrain)[_qp](c_comp, i) + (*_dW0p_broken_dstrain)[_qp](i, c_comp) ) * _grad_phi[_j][_qp](i) / _dt;
    }
  }

  // Contribution of the damage phase field to off diag Jacobian of temperature kernel
  if (jvar == _c_var)
    return 2.0 *  _plastic_factor * (1.0 - _c[_qp]) * ( _W0p[_qp] - _W0p_old[_qp] ) * _phi[_j][_qp] * _test[_i][_qp] / _dt;

  return - val * _test[_i][_qp];
}
