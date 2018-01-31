/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTETHERMALEXPANSIONEIGENSTRAINANISOTROPIC_H
#define COMPUTETHERMALEXPANSIONEIGENSTRAINANISOTROPIC_H

#include "ComputeThermalExpansionEigenstrainAnisotropicBase.h"
#include "DerivativeMaterialInterface.h"

class ComputeThermalExpansionEigenstrainAnisotropic;

template<>
InputParameters validParams<ComputeThermalExpansionEigenstrainAnisotropic>();

/**
 * ComputeThermalExpansionEigenstrain computes an eigenstrain for thermal expansion
 * with a constant expansion coefficient.
 */
class ComputeThermalExpansionEigenstrainAnisotropic : public ComputeThermalExpansionEigenstrainAnisotropicBase
{
public:
  ComputeThermalExpansionEigenstrainAnisotropic(const InputParameters & parameters);

protected:
    virtual void computeThermalStrain(std::vector<Real> & thermal_strain, std::vector<Real> & instantaneous_cte) override;

  const std::vector<Real> & _thermal_expansion_coeff;
};

#endif // COMPUTETHERMALEXPANSIONEIGENSTRAINANISOTROPIC_H
