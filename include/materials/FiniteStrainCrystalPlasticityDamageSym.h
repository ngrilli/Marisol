/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCRYSTALPLASTICITYDAMAGESYM_H
#define FINITESTRAINCRYSTALPLASTICITYDAMAGESYM_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticityDamageSym uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Lagrangian strain used to calculate damage (symmetric model)
 * Calculation of the broken plastic energy for temperature calculation
 */
class FiniteStrainCrystalPlasticityDamageSym;

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityDamageSym>();

class FiniteStrainCrystalPlasticityDamageSym : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCrystalPlasticityDamageSym(const InputParameters & parameters);

protected:
  /**
   * This function set variables for internal variable solve.
   */
  virtual void preSolveStatevar();

  /**
   * This function solves internal variables.
   */
  virtual void solveStatevar();

  /**
   * This function update internal variable after solve.
   */
  virtual void postSolveStatevar();


  // update slip system resistances and output slip increment
  virtual void updateGss();

  /**
   * Update elastic and plastic work
   */
  virtual void update_energies();

  /**
   * This function calculate stress residual.
   */
  virtual void calcResidual( RankTwoTensor & );

  virtual void getSlipIncrements();

  /**
   * This function calculate the elastic tangent moduli for preconditioner
   * xfac introduced for damage
   */
  virtual RankFourTensor elasticTangentModuli();

  const VariableValue & _c;

  /// Small number to avoid non-positive definiteness at or near complete damage
  Real _kdamage;

  // exponent in Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus
  const Real _Bulk_Modulus_Ref;

  // Von Neumann coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

  MaterialProperty<Real> & _W0e;
  MaterialProperty<Real> & _W0p;
  const MaterialProperty<Real> & _W0p_old;
  MaterialProperty<Real> & _W0p_broken;
  const MaterialProperty<Real> & _W0p_broken_old;
  MaterialProperty<RankTwoTensor> & _dstress_dc;
  MaterialProperty<RankTwoTensor> & _dW0e_dstrain;
  MaterialProperty<RankTwoTensor> & _dW0p_dstrain;
  MaterialProperty<RankTwoTensor> & _dW0p_broken_dstrain;
  MaterialProperty<RankTwoTensor> & _pk2_undamaged;
  MaterialProperty<RankTwoTensor> & _fe_out; // Elastic deformation gradient for output
  MaterialProperty<std::vector<Real>> & _slip_incr_out; // slip increment output

  Real _W0p_tmp;
  Real _W0p_tmp_old;
  Real _W0p_broken_tmp;
  Real _W0p_broken_tmp_old;

  std::vector<RankTwoTensor> _etens;
  std::vector<Real> _epos;
  std::vector<Real> _eigval;
  RankTwoTensor _eigvec;

};

#endif //FINITESTRAINCRYSTALPLASTICITYDAMAGESYM_H
