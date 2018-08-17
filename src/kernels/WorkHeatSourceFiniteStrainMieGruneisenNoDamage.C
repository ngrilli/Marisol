#include "WorkHeatSourceFiniteStrainMieGruneisenNoDamage.h"

template <>
InputParameters
validParams<WorkHeatSourceFiniteStrainMieGruneisenNoDamage>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Work heat source kernel"
                             "generic kernel for finite strain"
                             "sigma * dot(epsilon) term for energy conservation"
                             "Mie Gruneisen equation of state"
                             "(Menon, 2014) (Zhang, 2011)"
                             "no coupling with damage variable");
  params.addParam<std::string>("base_name",
                              "Optional parameter that allows the user to define "
                              "multiple mechanics material systems on the same "
                              "block, i.e. for multiple phases");
  params.addRequiredParam<Real>("G_Gruneisen", "Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS");
  params.addParam<MaterialPropertyName>(
      "specific_heat", "specific_heat", "Property name of the specific heat material property");
  params.addParam<MaterialPropertyName>(
      "density_name", "density", "Property name of the density material property");
  return params;
}

WorkHeatSourceFiniteStrainMieGruneisenNoDamage::WorkHeatSourceFiniteStrainMieGruneisenNoDamage(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _G_Gruneisen(getParam<Real>("G_Gruneisen")), // Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS
    _specific_heat(getMaterialProperty<Real>("specific_heat")),
    _density(getMaterialProperty<Real>("density_name")),
    _stress(getMaterialProperty<RankTwoTensor>(_base_name + "stress")), // Cauchy stress
    // In finite strain _strain_rate contains gradient calculated with the deformed mesh
    // therefore the mechanical work per unit volume is: _stress * _strain_rate = Cauchy_stress * dot(epsilon)
    _strain_rate(getMaterialProperty<RankTwoTensor>(_base_name + "strain_rate")), // strain rate = dot(epsilon)
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
WorkHeatSourceFiniteStrainMieGruneisenNoDamage::computeQpResidual()
{
  Real heat_source;

  heat_source = _stress[_qp].doubleContraction(_strain_rate[_qp]);

  return - heat_source * _test[_i][_qp];
}

Real
WorkHeatSourceFiniteStrainMieGruneisenNoDamage::computeQpJacobian()
{
  // approximate (small strain) Jacobian, first order in the strain rate
  // Jacobian = (d stress / dT) * strain_rate
  return _G_Gruneisen * _density[_qp] * _specific_heat[_qp] * _strain_rate[_qp].trace() * _phi[_j][_qp] * _test[_i][_qp];
}

Real
WorkHeatSourceFiniteStrainMieGruneisenNoDamage::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;

  val = 0.0;
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k]) {
      // approximate (small strain) off-diagonal Jacobian, first order in the strain rate
      // d/du_i ( stress * strain_rate ) = stress * d/du_i (strain_rate) = stress(j,h) * _grad_phi(j,h) * _test / dt
      for (unsigned int h = 0; h < _ndisp; ++h) {
        val += _stress[_qp](k,h) * _grad_phi[_j][_qp](h) * _test[_i][_qp] / _dt;
      }
    }
  }

  return val;
}
