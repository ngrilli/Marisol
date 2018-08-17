#ifndef THERMALEXPANSIONHEATSOURCEGASSOLID_H
#define THERMALEXPANSIONHEATSOURCEGASSOLID_H

#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class ThermalExpansionHeatSourceGasSolid;

template <>
InputParameters validParams<ThermalExpansionHeatSourceGasSolid>();

/**
 * This kernel calculates the heat source term corresponding to thermoelasticity
 * Birch Murnaghan equation of state for the solid phase
 * Mie Gruneisen equation of state for the gas phase
 * Luscher2017 finite strain formalism
 * isotropic elastic material
 * pressure average weighted with mass fraction
 */
class ThermalExpansionHeatSourceGasSolid : public HeatSource
{
public:
  ThermalExpansionHeatSourceGasSolid(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  // Mass fraction of the first specie
  const VariableValue & _mass_fraction_1;

  // Mass fraction of the second specie
  const VariableValue & _mass_fraction_2;

  // Murnaghan exponent in the equation of state of the solid phase
  const Real _n_Murnaghan;

  // Reference bulk modulus in the equation of state of the solid phase
  const Real _Bulk_Modulus_Ref_Solid;

  // Volumetric thermal expansion coefficient
  const Real _thermal_expansion;

  // Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS
  const Real _G_Gruneisen;

  // Cv in equation 7 in Menon 2014 to fit gas equation of state
  const Real _specific_heat_Gruneisen;

  // reference temperature with zero thermal expansion
  const Real _reference_temperature;

  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient; // deformation gradient

  // In finite strain _strain_rate contains gradient calculated with the deformed mesh
  // therefore the mechanical work per unit volume is: _stress * _strain_rate = Cauchy_stress * dot(epsilon)
  const MaterialProperty<RankTwoTensor> & _strain_rate; // strain rate dot(epsilon)

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // THERMALEXPANSIONHEATSOURCEGASSOLID_H
