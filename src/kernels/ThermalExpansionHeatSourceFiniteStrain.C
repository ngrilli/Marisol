#include "ThermalExpansionHeatSourceFiniteStrain.h"

template <>
InputParameters
validParams<ThermalExpansionHeatSourceFiniteStrain>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Thermal expansion heat source kernel"
                             "generic kernel for finite strain"
                             "Luscher2017 model and Amor damage model");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");
  params.addCoupledVar("c","Phase field damage variable: Used to indicate calculation of Off Diagonal Jacobian term");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  return params;
}

ThermalExpansionHeatSourceFiniteStrain::ThermalExpansionHeatSourceFiniteStrain(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _kdamage(getParam<Real>("kdamage")),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _thermal_expansion(getParam<Real>("thermal_expansion")), // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _deformation_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient
    _deformation_gradient_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")), // deformation gradient, previous timestep
    _fp(getMaterialPropertyByName<RankTwoTensor>(_base_name + "fp")), // Plastic deformation gradient
    _fp_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "fp")), // Plastic deformation gradient, previous timestep
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")), //elasticity tensor
    _c(coupledValue("c")),
    _c_coupled(isCoupled("c")),
    _c_var(_c_coupled ? coupled("c") : 0),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
ThermalExpansionHeatSourceFiniteStrain::computeQpResidual()
{
  Real heat_source, Je;
  Real Kb = _Bulk_Modulus_Ref;
  RankTwoTensor ce, ce_old, invce, iden, ee, ee_old, ee_rate, invce_ee_rate;
  RankTwoTensor fe, fe_old, inv_fp, inv_fp_old;

  // calculate elastic deformation gradient
  inv_fp = _fp[_qp].inverse();
  inv_fp_old = _fp_old[_qp].inverse();
  fe = _deformation_gradient[_qp] * inv_fp;
  fe_old = _deformation_gradient_old[_qp] * inv_fp_old;

  // preliminary calculation of deformation-related tensors
  Je = fe.det(); // Jacobian = relative volume
  ce = fe.transpose() * fe; // Ce = Cauchy tensor
  ce_old = fe_old.transpose() * fe_old; // Ce old = Cauchy tensor at previous time step
  invce = ce.inverse();
  iden.zero();
  iden.addIa(1.0);
  ee = 0.5 * (ce - iden);
  ee_old = 0.5 * (ce_old - iden);
  ee_rate = (ee - ee_old) / _dt;
  invce_ee_rate = invce * ee_rate;

  // Equation of state contribution (Psi_EOS in Luscher2017, equation 16)
  // heat source = - K alpha T exp(-n alpha (T - T0)) Je Tr(Ce^-1 dot(Ee))
  heat_source = - Kb * _thermal_expansion * _u[_qp]
              * std::exp(-_n_Murnaghan * _thermal_expansion * (_u[_qp] - _reference_temperature))
              * Je * invce_ee_rate.trace();

  if (Je > 1.0) {
    heat_source = heat_source * (1.0 - _c[_qp]) * (1.0 - _c[_qp]);
  }

  // Coupling contribution (Psi_cpl in Luscher2017, equation 15)
  // heat source = - alpha T / 3 exp(2/3 alpha (T-T0)) (dot(Ee) : C : I)
  //               + K alpha T Je^2/3 Tr(Ce^-1 dot(Ee)) exp(2/3 alpha (T-T0))
  RankTwoTensor thermal_coupling_tensor;
  thermal_coupling_tensor = _elasticity_tensor[_qp] * iden;
  thermal_coupling_tensor = ee_rate * thermal_coupling_tensor;
  heat_source += _thermal_expansion * _u[_qp]
              * std::exp((2.0/3.0) * _thermal_expansion * (_u[_qp] - _reference_temperature))
              * (Kb * std::pow(Je , 2.0/3.0) * invce_ee_rate.trace() - (1.0/3.0) * thermal_coupling_tensor.trace())
              * (1.0 - _c[_qp]) * (1.0 - _c[_qp]);

  return - heat_source * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceFiniteStrain::computeQpJacobian()
{
  Real Kb = _Bulk_Modulus_Ref;
  RankTwoTensor ce, ce_old, iden, ee, ee_old, ee_rate;
  RankTwoTensor fe, fe_old, inv_fp, inv_fp_old;

  // calculate elastic deformation gradient
  inv_fp = _fp[_qp].inverse();
  inv_fp_old = _fp_old[_qp].inverse();
  fe = _deformation_gradient[_qp] * inv_fp;
  fe_old = _deformation_gradient_old[_qp] * inv_fp_old;

  ce = fe.transpose() * fe; // Ce = Cauchy tensor
  ce_old = fe_old.transpose() * fe_old; // Ce old = Cauchy tensor at previous time step
  iden.zero();
  iden.addIa(1.0);
  ee = 0.5 * (ce - iden);
  ee_old = 0.5 * (ce_old - iden);
  ee_rate = (ee - ee_old) / _dt;

  // approximate (small strain) Jacobian, first order in the strain rate and temperature
  return (1.0/3.0) * Kb * _thermal_expansion * ee_rate.trace() * _phi[_j][_qp] * _test[_i][_qp];
}

Real
ThermalExpansionHeatSourceFiniteStrain::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  Real Kb = _Bulk_Modulus_Ref;
  RankTwoTensor ce, ce_old, iden, ee, ee_old, ee_rate;
  RankTwoTensor fe, fe_old, inv_fp, inv_fp_old;

  // calculate elastic deformation gradient
  inv_fp = _fp[_qp].inverse();
  inv_fp_old = _fp_old[_qp].inverse();
  fe = _deformation_gradient[_qp] * inv_fp;
  fe_old = _deformation_gradient_old[_qp] * inv_fp_old;

  ce = fe.transpose() * fe; // Ce = Cauchy tensor
  ce_old = fe_old.transpose() * fe_old; // Ce old = Cauchy tensor at previous time step
  iden.zero();
  iden.addIa(1.0);
  ee = 0.5 * (ce - iden);
  ee_old = 0.5 * (ce_old - iden);
  ee_rate = (ee - ee_old) / _dt;

  val = 0.0;
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k]) {
      // approximate (small strain) off-diagonal Jacobian, first order in the strain rate and temperature
      val = (1.0/3.0) * Kb * _thermal_expansion * _u[_qp] * _grad_phi[_j][_qp](k) * _test[_i][_qp] / _dt;
    }
  }

  if (jvar == _c_var) {
    val = - (2.0/3.0) * (1.0 - _c[_qp]) * Kb * _thermal_expansion * _u[_qp] * ee_rate.trace()
        * _phi[_j][_qp] * _test[_i][_qp];
  }

  return val;
}
