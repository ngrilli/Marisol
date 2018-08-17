#ifndef THERMALEXPANSIONHEATSOURCEBIRCHMURNAGHANISO_H
#define THERMALEXPANSIONHEATSOURCEBIRCHMURNAGHANISO_H

#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class ThermalExpansionHeatSourceBirchMurnaghanIso;

template <>
InputParameters validParams<ThermalExpansionHeatSourceBirchMurnaghanIso>();

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Luscher2017 model
 * Birch Murnaghan equation of state
 * isotropic elastic material
 */
class ThermalExpansionHeatSourceBirchMurnaghanIso : public HeatSource
{
public:
  ThermalExpansionHeatSourceBirchMurnaghanIso(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  // Murnaghan exponent in the equation of state
  const Real _n_Murnaghan;

  // Reference bulk modulus in the equation of state
  const Real _Bulk_Modulus_Ref;

  // Volumetric thermal expansion coefficient
  const Real _thermal_expansion;

  // reference temperature with zero thermal expansion
  const Real _reference_temperature;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient; // deformation gradient

  // In finite strain _strain_rate contains gradient calculated with the deformed mesh
  // therefore the mechanical work per unit volume is: _stress * _strain_rate = Cauchy_stress * dot(epsilon)
  const MaterialProperty<RankTwoTensor> & _strain_rate; // strain rate dot(epsilon)

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // THERMALEXPANSIONHEATSOURCEBIRCHMURNAGHANISO_H
