#ifndef THERMALEXPANSIONHEATSOURCEFINITESTRAIN_H
#define THERMALEXPANSIONHEATSOURCEFINITESTRAIN_H

#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class ThermalExpansionHeatSourceFiniteStrain;

template <>
InputParameters validParams<ThermalExpansionHeatSourceFiniteStrain>();

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Luscher2017 model and Amor damage model
 */
class ThermalExpansionHeatSourceFiniteStrain : public HeatSource
{
public:
  ThermalExpansionHeatSourceFiniteStrain(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  // residual stiffness at c = 1
  const Real _kdamage;

  // Murnaghan exponent in the equation of state
  const Real _n_Murnaghan;

  // Reference bulk modulus in the equation of state
  const Real _Bulk_Modulus_Ref;

  // Volumetric thermal expansion coefficient
  const Real _thermal_expansion;

  // reference temperature with zero thermal expansion
  const Real _reference_temperature;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient; // deformation gradient
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old; // deformation gradient, previous timestep

  const MaterialProperty<RankTwoTensor> & _fp; // plastic deformation gradient
  const MaterialProperty<RankTwoTensor> & _fp_old; // plastic deformation gradient, previous timestep

  const MaterialProperty<RankFourTensor> & _elasticity_tensor; //elasticity tensor

  // damage phase field
  const VariableValue & _c;
  const bool _c_coupled;
  const unsigned int _c_var;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // THERMALEXPANSIONHEATSOURCEFINITESTRAIN_H
