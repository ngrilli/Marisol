//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatConductionMaterial2species.h"
#include "Function.h"

#include "libmesh/quadrature.h"

registerMooseObject("HeatConductionApp", HeatConductionMaterial2species);

template <>
InputParameters
validParams<HeatConductionMaterial2species>()
{
  InputParameters params = validParams<Material>();

  params.addCoupledVar("temp", "Coupled Temperature");
  params.addCoupledVar("mass_fraction_1", "Mass fraction of the first specie");
  params.addCoupledVar("mass_fraction_2", "Mass fraction of the second specie");
// third specie is (1.0 - mass_fraction_1 - mass_fraction_2)

  params.addParam<Real>("thermal_conductivity", "The thermal conductivity value");
  params.addParam<FunctionName>("thermal_conductivity_temperature_function",
                                "",
                                "Thermal conductivity as a function of temperature.");

  params.addParam<Real>("specific_heat_1", "The specific heat value of the first specie");
  params.addParam<FunctionName>(
      "specific_heat_temperature_function_1", "", "Specific heat of the first specie as a function of temperature.");
  params.addParam<Real>("specific_heat_2", "The specific heat value of the second specie");
  params.addParam<FunctionName>(
      "specific_heat_temperature_function_2", "", "Specific heat of the second specie as a function of temperature.");
  params.addParam<Real>("specific_heat_3", "The specific heat value of the third specie");
  params.addParam<FunctionName>(
      "specific_heat_temperature_function_3", "", "Specific heat of the third specie as a function of temperature.");

  params.addClassDescription("General-purpose material model for heat conduction: 2 species");

  return params;
}

HeatConductionMaterial2species::HeatConductionMaterial2species(const InputParameters & parameters)
  : Material(parameters),

    _has_temp(isCoupled("temp")),
    _temperature(_has_temp ? coupledValue("temp") : _zero),
    _has_mass_fraction_1(isCoupled("mass_fraction_1")),
    _mass_fraction_1(_has_mass_fraction_1 ? coupledValue("mass_fraction_1") : _zero),
    _has_mass_fraction_2(isCoupled("mass_fraction_2")),
    _mass_fraction_2(_has_mass_fraction_2 ? coupledValue("mass_fraction_2") : _zero),
    // mass_fraction_3 = 1.0 - mass_fraction_1 - mass_fraction_2

    _my_thermal_conductivity(
        isParamValid("thermal_conductivity") ? getParam<Real>("thermal_conductivity") : 0),
    _my_specific_heat_1(isParamValid("specific_heat_1") ? getParam<Real>("specific_heat_1") : 0),
    _my_specific_heat_2(isParamValid("specific_heat_2") ? getParam<Real>("specific_heat_2") : 0),
    _my_specific_heat_3(isParamValid("specific_heat_3") ? getParam<Real>("specific_heat_3") : 0),

    _thermal_conductivity(declareProperty<Real>("thermal_conductivity")),
    _thermal_conductivity_dT(declareProperty<Real>("thermal_conductivity_dT")),
    _thermal_conductivity_temperature_function(
        getParam<FunctionName>("thermal_conductivity_temperature_function") != ""
            ? &getFunction("thermal_conductivity_temperature_function")
            : NULL),

    _specific_heat(declareProperty<Real>("specific_heat")),
    _specific_heat_temperature_function_1(
        getParam<FunctionName>("specific_heat_temperature_function_1") != ""
            ? &getFunction("specific_heat_temperature_function_1")
            : NULL),
    _specific_heat_temperature_function_2(
        getParam<FunctionName>("specific_heat_temperature_function_2") != ""
            ? &getFunction("specific_heat_temperature_function_2")
            : NULL),
    _specific_heat_temperature_function_3(
        getParam<FunctionName>("specific_heat_temperature_function_3") != ""
            ? &getFunction("specific_heat_temperature_function_3")
            : NULL)
{
  if (_thermal_conductivity_temperature_function && !_has_temp)
  {
    mooseError("Must couple with temperature if using thermal conductivity function");
  }
  if (isParamValid("thermal_conductivity") && _thermal_conductivity_temperature_function)
  {
    mooseError(
        "Cannot define both thermal conductivity and thermal conductivity temperature function");
  }
  if (_specific_heat_temperature_function_1 && !_has_temp)
  {
    mooseError("Must couple with temperature if using specific heat function");
  }
  if (isParamValid("specific_heat_1") && _specific_heat_temperature_function_1)
  {
    mooseError("Cannot define both specific heat and specific heat temperature function");
  }
}

void
HeatConductionMaterial2species::computeProperties()
{
  Real mass_fraction_1 = std::min(_mass_fraction_1[_qp] , 1.0);
  mass_fraction_1 = std::max(mass_fraction_1 , 0.0);

  Real mass_fraction_2 = std::min(_mass_fraction_2[_qp] , 1.0);
  mass_fraction_2 = std::max(mass_fraction_2 , 0.0);

  // sum of mole fractions is one
  Real mass_fraction_3 = 1.0 - mass_fraction_1 - mass_fraction_2;

  for (unsigned int qp(0); qp < _qrule->n_points(); ++qp)
  {
    Real qp_temperature = 0;
    if (_has_temp)
    {
      qp_temperature = _temperature[qp];
      if (_temperature[qp] < 0)
      {
        std::stringstream msg;
        msg << "WARNING:  In HeatConductionMaterial:  negative temperature!\n"
            << "\tResetting to zero.\n"
            << "\t_qp: " << qp << "\n"
            << "\ttemp: " << _temperature[qp] << "\n"
            << "\telem: " << _current_elem->id() << "\n"
            << "\tproc: " << processor_id() << "\n";
        mooseWarning(msg.str());
        qp_temperature = 0;
      }
    }
    if (_thermal_conductivity_temperature_function)
    {
      Point p;
      _thermal_conductivity[qp] =
          _thermal_conductivity_temperature_function->value(qp_temperature, p);
      _thermal_conductivity_dT[qp] =
          _thermal_conductivity_temperature_function->timeDerivative(qp_temperature, p);
    }
    else
    {
      _thermal_conductivity[qp] = _my_thermal_conductivity;
      _thermal_conductivity_dT[qp] = 0;
    }

    if (_specific_heat_temperature_function_1) // mixture not implemented for function input
    {
      Point p;
      _specific_heat[qp] = _specific_heat_temperature_function_1->value(qp_temperature, p);
    }
    else
    {
      _specific_heat[qp] = _my_specific_heat_1 * mass_fraction_1
                         + _my_specific_heat_2 * mass_fraction_2
                         + _my_specific_heat_3 * mass_fraction_3;
    }
  }
}
