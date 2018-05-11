//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef HEATCONDUCTIONMATERIAL2SPECIES_H
#define HEATCONDUCTIONMATERIAL2SPECIES_H

#include "Material.h"

// Forward Declarations
class HeatConductionMaterial2species;
class Function;

template <>
InputParameters validParams<HeatConductionMaterial2species>();

/**
 * General-purpose material model for heat conduction: 2 species
 */
class HeatConductionMaterial2species : public Material
{
public:
  HeatConductionMaterial2species(const InputParameters & parameters);

protected:
  virtual void computeProperties();

  const bool _has_temp;
  const VariableValue & _temperature;

  const bool _has_mole_fraction_1;
  const VariableValue & _mole_fraction_1;

  const bool _has_mole_fraction_2;
  const VariableValue & _mole_fraction_2;

  const Real _molar_mass_1;
  const Real _molar_mass_2;
  const Real _molar_mass_3;

  const Real _my_thermal_conductivity;
  const Real _my_specific_heat_1;
  const Real _my_specific_heat_2;
  const Real _my_specific_heat_3;

  MaterialProperty<Real> & _thermal_conductivity;
  MaterialProperty<Real> & _thermal_conductivity_dT;
  Function * _thermal_conductivity_temperature_function;

  MaterialProperty<Real> & _specific_heat;
  Function * _specific_heat_temperature_function_1;
  Function * _specific_heat_temperature_function_2;
  Function * _specific_heat_temperature_function_3;
};

#endif // HEATCONDUCTIONMATERIAL2SPECIES_H
