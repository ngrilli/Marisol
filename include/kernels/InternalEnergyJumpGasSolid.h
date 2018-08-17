#ifndef INTERNALENERGYJUMPGASSOLID_H
#define INTERNALENERGYJUMPGASSOLID_H

#include "DerivativeMaterialInterface.h"
#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class InternalEnergyJumpGasSolid;

template <>
InputParameters validParams<InternalEnergyJumpGasSolid>();

/**
 * Heat source kernel
 * for the solid-gas mixture
 * to take into account the term
 * (e_{BM} - e_{MG}) dot(Y_3)
 * representing the jump in internal energy
 * when solid and intermediates transform into gas
 * must be added to the term ThermalExpansionHeatSourceGasSolid
 * to be thermodynamically consistent
 * Birch Murnaghan equation of state for the solid phase
 * Mie Gruneisen equation of state for the gas phase
 * Luscher2017 finite strain formalism
 * isotropic elastic material
 * pressure average weighted with mass fraction
 */
class InternalEnergyJumpGasSolid : public DerivativeMaterialInterface<HeatSource>
{
public:
  InternalEnergyJumpGasSolid(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  // Mass fraction of the second specie (intermediates)
  const VariableValue & _mass_fraction_2;

  /// MOOSE variable number for the mass fraction variable
  // of the second specie (intermediates)
  const unsigned int _mass_fraction_2_var;

  /// name of mass fraction 2 variable (intermediates)
  VariableName _mass_fraction_2_name;

  // Murnaghan exponent in the equation of state of the solid phase
  const Real _n_Murnaghan;

  // Reference bulk modulus in the equation of state of the solid phase
  const Real _Bulk_Modulus_Ref_Solid;

  // KT0 prime correction (Menikoff 2001, Yoo 1998)
  const Real _Bulk_Modulus_Cor_Solid;

  // Volumetric thermal expansion coefficient
  const Real _thermal_expansion;

  // Gruneisen coefficient G (or Gamma) in Mie-Gruneisen EOS
  const Real _G_Gruneisen;

  // Us-Up slope in Mie-Gruneisen EOS
  const Real _s_UsUp;

  // reference bulk modulus = rho_0 c_0^2 (equation 12 in Menon 2014)
  const Real _bulk_modulus_ref_gas;

  // Cv in equation 7 in Menon 2014 to fit gas equation of state
  const Real _specific_heat_Gruneisen;

  // reference temperature with zero thermal expansion
  const Real _reference_temperature;

  const MaterialProperty<Real> & _density;

  // _mass_fraction_rate is dot(Y_3) from arrhenius_material_Y2_omega2 (positive prefactor)
  const MaterialProperty<Real> & _mass_fraction_rate;
  const MaterialProperty<Real> & _dmass_fraction_rate_dtemperature;
  const MaterialProperty<Real> & _dmass_fraction_rate_dmass_fraction;

  const MaterialProperty<RankTwoTensor> & _deformation_gradient; // deformation gradient

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // INTERNALENERGYJUMPGASSOLID_H
