#ifndef PLASTICHEATINGSOURCEDAMAGE_H
#define PLASTICHEATINGSOURCEDAMAGE_H

#include "HeatSource.h"
#include "RankTwoTensor.h"

// Forward Declarations
class PlasticHeatingSourceDamage;

template <>
InputParameters validParams<PlasticHeatingSourceDamage>();

/**
 * This kernel calculates the heat source term corresponding to plastic deformation
 * plasticity with damage (Miehe 2016)
 */
class PlasticHeatingSourceDamage : public HeatSource
{
public:
  PlasticHeatingSourceDamage(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  /// Pre-Factor of the plastic energy contribution to damage evolution
  Real _plastic_factor;

  /// Contribution of plastic strain energy to damage evolution
  const MaterialProperty<Real> & _W0p; // plastic energy (unbroken)
  const MaterialProperty<Real> & _W0p_old; // old plastic energy (unbroken)
  const MaterialProperty<Real> & _W0p_broken; // plastic energy (broken)
  const MaterialProperty<Real> & _W0p_broken_old; // old plastic energy (broken)

  /// Variation of plastic strain energy driving damage evolution with strain
  const MaterialProperty<RankTwoTensor> * _dW0p_dstrain;
  const MaterialProperty<RankTwoTensor> * _dW0p_broken_dstrain;

  /// order parameter for damage
  const VariableValue & _c;
  const unsigned int _c_var;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // PLASTICHEATINGSOURCEDAMAGE_H
