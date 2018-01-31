/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputeThermalExpansionEigenstrainAnisotropicBase.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"
#include "ElementPropertyReadFile.h"

template<>
InputParameters validParams<ComputeThermalExpansionEigenstrainAnisotropicBase>()
{
  InputParameters params = validParams<ComputeEigenstrainBase>();
  params.addCoupledVar("temperature", "Coupled temperature");
  params.addDeprecatedParam<Real>("stress_free_reference_temperature", "Reference temperature for thermal eigenstrain calculation", "'stress_free_temperature' has replaced this parameter");
  params.addParam<Real>("stress_free_temperature", "Reference temperature for thermal eigenstrain calculation");
  params.addParam<UserObjectName>("read_prop_user_object_eigenstrain","The ElementReadPropertyFile GeneralUserObject to read element specific property values from file");
    
  return params;
}

ComputeThermalExpansionEigenstrainAnisotropicBase::ComputeThermalExpansionEigenstrainAnisotropicBase(const InputParameters & parameters) :
    DerivativeMaterialInterface<ComputeEigenstrainBase>(parameters),
    _temperature(coupledValue("temperature")),
    _deigenstrain_dT(declarePropertyDerivative<RankTwoTensor>(_eigenstrain_name, getVar("temperature",0)->name())),
   _read_prop_user_object_eigenstrain(isParamValid("read_prop_user_object_eigenstrain") ? & getUserObject<ElementPropertyReadFile>("read_prop_user_object_eigenstrain") : NULL),
   _Euler_angles_mat_prop(declareProperty<RealVectorValue>("Euler_angles_eigenstrain")),
   _crysrot_eigenstrain(declareProperty<RankTwoTensor>("crysrot_eigenstrain")),
   _R(_Euler_angles_eigenstrain)

{
  if (isParamValid("stress_free_temperature"))
    _stress_free_temperature = getParam<Real>("stress_free_temperature");
  else if (isParamValid("stress_free_reference_temperature"))
    _stress_free_temperature = getParam<Real>("stress_free_reference_temperature");
  else
    mooseError("Please specify 'stress_free_temperature'.");
}

void
ComputeThermalExpansionEigenstrainAnisotropicBase::assignEulerAngles()
{
    if (_read_prop_user_object_eigenstrain)
    {
        _Euler_angles_mat_prop[_qp](0) = _read_prop_user_object_eigenstrain->getData(_current_elem, 0);
        _Euler_angles_mat_prop[_qp](1) = _read_prop_user_object_eigenstrain->getData(_current_elem, 1);
        _Euler_angles_mat_prop[_qp](2) = _read_prop_user_object_eigenstrain->getData(_current_elem, 2);
      }
    else
        _Euler_angles_mat_prop[_qp] = _Euler_angles_eigenstrain;
}


void
ComputeThermalExpansionEigenstrainAnisotropicBase::computeQpEigenstrain()
{
    std::vector<Real> thermal_strain;
    thermal_strain.resize(3);
    std::vector<Real> instantaneous_cte;
    instantaneous_cte.resize(3);

  computeThermalStrain(thermal_strain, instantaneous_cte);
    //Properties assigned at the beginning of every call to material calculation
    
    assignEulerAngles();
    
    _R.update(_Euler_angles_mat_prop[_qp]);
    
  _eigenstrain[_qp].zero();
  _eigenstrain[_qp](0,0) = thermal_strain[0];
  _eigenstrain[_qp](1,1) = thermal_strain[1];
  _eigenstrain[_qp](2,2) = thermal_strain[2];
  _crysrot_eigenstrain[_qp] = _R.transpose();
  _eigenstrain[_qp].rotate(_crysrot_eigenstrain[_qp]);
    
  _deigenstrain_dT[_qp].zero();
  _deigenstrain_dT[_qp](0,0) = instantaneous_cte[0];
  _deigenstrain_dT[_qp](1,1) = instantaneous_cte[1];
  _deigenstrain_dT[_qp](2,2) = instantaneous_cte[2];
    
  _deigenstrain_dT[_qp].rotate(_crysrot_eigenstrain[_qp]);

}

