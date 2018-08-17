#ifndef WORKHEATSOURCEFINITESTRAINMIEGRUNEISENNODAMAGE_H
#define WORKHEATSOURCEFINITESTRAINMIEGRUNEISENNODAMAGE_H

#include "HeatSource.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class WorkHeatSourceFiniteStrainMieGruneisenNoDamage;

template <>
InputParameters validParams<WorkHeatSourceFiniteStrainMieGruneisenNoDamage>();

/**
 * This kernel calculates the heat source term corresponding to
   mechanical work:
   generic kernel for finite strain
   sigma * dot(epsilon) term for energy conservation
   Mie Gruneisen equation of state
   (Menon, 2014) (Zhang, 2011)
   no coupling with damage variable
 */
class WorkHeatSourceFiniteStrainMieGruneisenNoDamage : public HeatSource
{
public:
  WorkHeatSourceFiniteStrainMieGruneisenNoDamage(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  // Gruneisen G (or Gamma) parameter (Menon, 2014)
  const Real _G_Gruneisen;

  const MaterialProperty<Real> & _specific_heat;
  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _stress; // Cauchy stress

  // In finite strain _strain_rate contains gradient calculated with the deformed mesh
  // therefore the mechanical work per unit volume is: _stress * _strain_rate = Cauchy_stress * dot(epsilon)
  const MaterialProperty<RankTwoTensor> & _strain_rate; // strain rate dot(epsilon)

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // WORKHEATSOURCEFINITESTRAINMIEGRUNEISENNODAMAGE_H
