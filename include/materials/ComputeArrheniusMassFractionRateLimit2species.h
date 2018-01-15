/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEARRHENIUSMASSFRACTIONRATELIMIT2SPECIES_H
#define COMPUTEARRHENIUSMASSFRACTIONRATELIMIT2SPECIES_H

#include "DerivativeMaterialInterface.h"
#include "Material.h"

/**
 * ComputeArrheniusMassFractionRateLimit2species computes the mass fraction rate
 * for 2 species given by chemical reactions, upper threshold rate is used
 * specie 1 is the kernel variable
 * species 2 allows to calculate the Off diagonal Jacobian
 */
class ComputeArrheniusMassFractionRateLimit2species : public DerivativeMaterialInterface<Material>
{
public:
  ComputeArrheniusMassFractionRateLimit2species(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// optional parameter that allows multiple mechanics materials to be defined
  std::string _base_name;

  /// Auxiliary variable: mass_fraction_1 (kernel variable)
  const VariableValue & _mass_fraction_1;

  /// MOOSE variable number for the mass fraction 1 variable
  const unsigned int _mass_fraction_1_var;

  /// name of mass fraction 1 variable
  VariableName _mass_fraction_1_name;

  /// Auxiliary variable: mass_fraction_2
  const VariableValue & _mass_fraction_2;

  /// MOOSE variable number for the mass fraction 2 variable
  const unsigned int _mass_fraction_2_var;

  /// name of mass fraction 2 variable
  VariableName _mass_fraction_2_name;

  /// Auxiliary variable: temperature
  const VariableValue & _temperature;

  /// MOOSE variable number for the temperature variable
  const unsigned int _temperature_var;

  /// name of temperature variable
  VariableName _temperature_name;

  Real _exponential_prefactor;
  Real _exponential_coefficient;
  Real _exponential_factor;

  /// Upper threshold for the reaction rate
  Real _rate_limit;

  /// Exponential factor of the mass fraction of the species
  Real _nu_1;
  Real _nu_2;

  /// computed property: mass fraction rate
  MaterialProperty<Real> & _mass_fraction_rate;

  /// d(mass_fraction_rate)/d(T) for Off diagonal Jacobian
  MaterialProperty<Real> & _dmass_fraction_rate_dtemperature;

  /// d(mass_fraction_rate)/d(mass_fraction_1) for diagonal Jacobian
  MaterialProperty<Real> & _dmass_fraction_rate_dmass_fraction_1;

  /// d(mass_fraction_rate)/d(mass_fraction_2) for Off diagonal Jacobian
  MaterialProperty<Real> & _dmass_fraction_rate_dmass_fraction_2;
};

#endif // COMPUTEARRHENIUSMASSFRACTIONRATELIMIT2SPECIES_H
