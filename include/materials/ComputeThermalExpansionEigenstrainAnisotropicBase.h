/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTETHERMALEXPANSIONEIGENSTRAINANISOTROPICBASE_H
#define COMPUTETHERMALEXPANSIONEIGENSTRAINANISOTROPICBASE_H

#include "ComputeEigenstrainBase.h"
#include "DerivativeMaterialInterface.h"

#include "ElementPropertyReadFile.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"

class ComputeThermalExpansionEigenstrainAnisotropicBase;
class RankTwoTensor;

template<>
InputParameters validParams<ComputeThermalExpansionEigenstrainAnisotropicBase>();

/**
 * ComputeThermalExpansionEigenstrainBase is a base class for all models that
 * compute eigenstrains due to thermal expansion of a material.
 */
class ComputeThermalExpansionEigenstrainAnisotropicBase : public DerivativeMaterialInterface<ComputeEigenstrainBase>
{
public:
  ComputeThermalExpansionEigenstrainAnisotropicBase(const InputParameters & parameters);

protected:
  virtual void computeQpEigenstrain() override;
  virtual void computeThermalStrain(std::vector<Real> & thermal_strain, std::vector<Real> & instantaneous_cte) = 0;
  virtual void assignEulerAngles();

  const VariableValue & _temperature;
  MaterialProperty<RankTwoTensor> & _deigenstrain_dT;
  Real _stress_free_temperature;
  const ElementPropertyReadFile * _read_prop_user_object_eigenstrain;
    
  MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;
    
    /// Crystal Rotation Matrix
  MaterialProperty<RankTwoTensor> & _crysrot_eigenstrain;
    
    /// Rotation matrix
  RotationTensor _R;
 RealVectorValue _Euler_angles_eigenstrain;

};

#endif // COMPUTETHERMALEXPANSIONEIGENSTRAINANISOTROPICBASE_H
