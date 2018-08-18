/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainCrystalPlasticityDamageAmorEOS.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityDamageAmorEOS>()
{
  InputParameters params = validParams<FiniteStrainCrystalPlasticity>();
  params.addClassDescription("Crystal Plasticity base class. Damage. Amor 2009."
                             "calculation of the broken plastic energy for temperature calculation"
                             "third-order Birch Murnaghan equation of state"
                             "finite strain formalism as in Luscher2017");
  params.addRequiredCoupledVar("c","Order parameter for damage");
  params.addRequiredCoupledVar("temp","Temperature");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("bulk_modulus_cor", "KT0 prime correction (Menikoff 2001)");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");

  return params;
}

FiniteStrainCrystalPlasticityDamageAmorEOS::FiniteStrainCrystalPlasticityDamageAmorEOS(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
    _c(coupledValue("c")),
    _temp(coupledValue("temp")),
    _kdamage(getParam<Real>("kdamage")),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _Bulk_Modulus_Cor(getParam<Real>("bulk_modulus_cor")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _thermal_expansion(getParam<Real>("thermal_expansion")), // volumetric thermal exapnsion coeffcient, as in Austin Barton 2015
    _reference_temperature(getParam<Real>("reference_temperature")), // reference temperature, as in Luscher2017
    _W0e(declareProperty<Real>("W0e")), // elastic energy
    _W0p(declareProperty<Real>("W0p")), // plastic energy (unbroken)
    _W0p_old(getMaterialPropertyOld<Real>("W0p")), // plastic energy of previous increment (unbroken)
    _W0p_broken(declareProperty<Real>("W0p_broken")), // plastic energy (broken)
    _W0p_broken_old(getMaterialPropertyOld<Real>("W0p_broken")), // plastic energy of previous increment (broken)
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _dW0e_dstrain(declareProperty<RankTwoTensor>("dW0e_dstrain")),
    _dW0p_dstrain(declareProperty<RankTwoTensor>("dW0p_dstrain")),
    _dW0p_broken_dstrain(declareProperty<RankTwoTensor>("dW0p_broken_dstrain")), // for the Jacobian calculation in the temperature kernel
    _pk2_undamaged(declareProperty<RankTwoTensor>("pk2_undamaged")), // undamaged 2nd Piola Kirchoff Stress
    _fe_out(declareProperty<RankTwoTensor>("fe_out")), // Elastic deformation gradient for output
    _slip_incr_out(declareProperty<std::vector<Real>>("slip_incr_out")) // slip system strain increment for output
{
}

void
FiniteStrainCrystalPlasticityDamageAmorEOS::preSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss_tmp = _gss_old[_qp];
    _W0p_tmp = _W0p_old[_qp];
    _W0p_broken_tmp = _W0p_broken_old[_qp];
    _accslip_tmp_old = _acc_slip_old[_qp];
  }
  else
  {
    if (_first_step_iter)
    {
      _gss_tmp = _gss_tmp_old = _gss_old[_qp];
      _W0p_tmp = _W0p_tmp_old = _W0p_old[_qp];
      _W0p_broken_tmp = _W0p_broken_tmp_old = _W0p_broken_old[_qp];
      _accslip_tmp_old = _acc_slip_old[_qp];
    }
    else
    {
      _gss_tmp = _gss_tmp_old;
      _W0p_tmp = _W0p_tmp_old;
      _W0p_broken_tmp = _W0p_broken_tmp_old;
    }
  }
}

void
FiniteStrainCrystalPlasticityDamageAmorEOS::solveStatevar()
{
  Real gmax, gdiff;
  unsigned int iterg;
  std::vector<Real> gss_prev(_nss);

  gmax = 1.1 * _gtol;
  iterg = 0;

  while (gmax > _gtol && iterg < _maxiterg) // Check for slip system resistance update tolerance
  {
    preSolveStress();
    solveStress();
    if (_err_tol)
      return;

    update_energies(); // update elastic and plastic energies

    postSolveStress(); // Update _fp_old_inv = _fp_old

    gss_prev = _gss_tmp;

    update_slip_system_resistance(); // Update slip system resistance

    gmax = 0.0;
    for (unsigned i = 0; i < _nss; ++i)
    {
      gdiff = std::abs(gss_prev[i] - _gss_tmp[i]); // Calculate increment size

      if (gdiff > gmax)
        gmax = gdiff;
    }
    iterg++;
  }

  if (iterg == _maxiterg)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainCrystalPLasticity: Hardness Integration error gmax", gmax, "\n");
#endif
    _err_tol = true;
  }
}

void
FiniteStrainCrystalPlasticityDamageAmorEOS::postSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss[_qp] = _gss_tmp;
    _W0p[_qp] = _W0p_tmp;
    _W0p_broken[_qp] = _W0p_broken_tmp;
    _acc_slip[_qp] = _accslip_tmp;
  }
  else
  {
    if (_last_step_iter)
    {
      _gss[_qp] = _gss_tmp;
      _W0p[_qp] = _W0p_tmp;
      _W0p_broken[_qp] = _W0p_broken_tmp;
      _acc_slip[_qp] = _accslip_tmp;
    }
    else
    {
      _gss_tmp_old = _gss_tmp;
      _W0p_tmp_old = _W0p_tmp;
      _W0p_broken_tmp_old = _W0p_broken_tmp;
      _accslip_tmp_old = _accslip_tmp;
    }
  }
}

/**
 * Old function to update slip system resistances.
 * Kept to avoid code break at computeQpstress
 * output slip increment
 */
void
FiniteStrainCrystalPlasticityDamageAmorEOS::updateGss()
{
  DenseVector<Real> hb(_nss);
  Real qab;

  Real a = _hprops[4]; // Kalidindi

  _slip_incr_out[_qp].resize(_nss);

  _accslip_tmp = _accslip_tmp_old;
  for (unsigned int i = 0; i < _nss; ++i)
    _accslip_tmp += std::abs(_slip_incr(i));

  for (unsigned int i = 0; i < _nss; ++i)
    _slip_incr_out[_qp][i] = _slip_incr(i);

  // Real val = std::cosh(_h0 * _accslip_tmp / (_tau_sat - _tau_init)); // Karthik
  // val = _h0 * std::pow(1.0/val,2.0); // Kalidindi

  for (unsigned int i = 0; i < _nss; ++i)
    // hb(i)=val;
    hb(i) = _h0 * std::pow(std::abs(1.0 - _gss_tmp[i] / _tau_sat), a) *
            copysign(1.0, 1.0 - _gss_tmp[i] / _tau_sat);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    if (_max_substep_iter == 1) // No substepping
      _gss_tmp[i] = _gss_old[_qp][i];
    else
      _gss_tmp[i] = _gss_tmp_old[i];

    for (unsigned int j = 0; j < _nss; ++j)
    {
      unsigned int iplane, jplane;
      iplane = i / 3;
      jplane = j / 3;

      if (iplane == jplane) // Kalidindi
        qab = 1.0;
      else
        qab = _r;

      _gss_tmp[i] += qab * hb(j) * std::abs(_slip_incr(j));
      _dgss_dsliprate(i, j) = qab * hb(j) * copysign(1.0, _slip_incr(j)) * _dt;
    }
  }
}

// Update slip system resistance and elastic and plastic work
void
FiniteStrainCrystalPlasticityDamageAmorEOS::update_energies()
{
  RankTwoTensor cauchy_stress_undamaged, cauchy_stress, WpToTrace, WpBrokenToTrace, invFe;
  Real detFe;

  if (_max_substep_iter == 1) //No substepping
  {
    _W0p_tmp = _W0p_old[_qp];
    _W0p_broken_tmp = _W0p_broken_old[_qp];
  }
  else
  {
    _W0p_tmp = _W0p_tmp_old;
    _W0p_broken_tmp = _W0p_broken_tmp_old;
  }

  // Update elastic and plastic work
  detFe = _fe.det();
  invFe = _fe.inverse();

  // _pk2[_qp] is the updated piola-kirchoff
  cauchy_stress = _fe * _pk2[_qp] * _fe.transpose()/detFe;
  cauchy_stress_undamaged = _fe * _pk2_undamaged[_qp] * _fe.transpose()/detFe;

  WpBrokenToTrace = cauchy_stress * _fe * ( _fp_inv.inverse() - _fp_old_inv.inverse() ) * _fp_inv * invFe * detFe;
  WpToTrace = cauchy_stress_undamaged * _fe * ( _fp_inv.inverse() - _fp_old_inv.inverse() ) * _fp_inv * invFe * detFe;

  _W0p_broken_tmp += WpBrokenToTrace.trace();
  _W0p_tmp += WpToTrace.trace();

  _dW0p_broken_dstrain[_qp].zero();
  _dW0p_dstrain[_qp].zero();
  //_dstress_dc[_qp] = -cauchy_stress_undamaged * (2.0 * (1.0 - c));

}

void
FiniteStrainCrystalPlasticityDamageAmorEOS::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, invce, ee, ce_pk2, eqv_slip_incr, pk2_new;
  Real trD, Kb, KT0prime, Je;
  Real c = _c[_qp];
  Real temp = _temp[_qp];
  Real xfac = Utility::pow<2>(1.0-c) + _kdamage;

  // Anisotropic elasticity (Luscher2017)
  // Kb = K in Luscher2017
  // Kb = (1/9) I : C : I
  //Real Kb = 0.0;

  //for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //  for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  //    Kb +=  _elasticity_tensor[_qp](i, i, j, j);

  //Kb = (1.0 / 9.0) * Kb;

  // or reference bulk modulus can be an input
  Kb = _Bulk_Modulus_Ref;
  KT0prime = _Bulk_Modulus_Cor;

  iden.zero();
  iden.addIa(1.0);

  _fe = _dfgrd_tmp * _fp_prev_inv; // _fp_inv  ==> _fp_prev_inv

  ce = _fe.transpose() * _fe;
  ce_pk2 = ce * _pk2_tmp;
  ce_pk2 = ce_pk2 / _fe.det();

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _nss; ++i)
    _tau(i) = ce_pk2.doubleContraction(_s0[i]);

  getSlipIncrements(); // Calculate dslip,dslipdtau

  if (_err_tol)
    return;

  eqv_slip_incr.zero();
  for (unsigned int i = 0; i < _nss; ++i)
    eqv_slip_incr += _s0[i] * _slip_incr(i);

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;
  _fe_out[_qp] = _fe; // Elastic deformation gradient for output

  ce = _fe.transpose() * _fe;
  invce = ce.inverse();
  ee = 0.5 * ( ce - iden );
  Je = _fe.det(); // Jacobian = relative volume

  Real Peos, Peos_pos, Peos_neg, V0V;

  // Cauchy pressure, third-order Birch Murnaghan (Menikoff2001, Yoo1998)
  V0V = 1.0 / Je; // relative volume
  Peos = 1.5 * Kb * (std::pow(V0V , 7.0/3.0) - std::pow(V0V , 5.0/3.0));
  Peos = Peos * (1.0 + 0.75 * (KT0prime - 4.0) * (std::pow(V0V , 2.0/3.0) - 1.0));
  Peos = - Peos; // negative stress in compression
  Peos += (Kb / _n_Murnaghan) * (std::exp(-_n_Murnaghan * _thermal_expansion * (temp - _reference_temperature)) - 1.0);

  Peos_pos = (std::abs(Peos) + Peos) / 2.0; // volumetric expansion
  Peos_neg = (std::abs(Peos) - Peos) / 2.0; // volumetric compression

  // positive and negative volumetric stresses (equation 17 in Luscher2017) + damage
  pk2_new = Je * (Peos_pos * xfac - Peos_neg) * invce;

  // Undamaged second piola-kirchoff stress to calculate undamaged plastic work
  _pk2_undamaged[_qp] = Je * Peos * invce;

  RankTwoTensor thermal_eigenstrain;
  // thermal eigenstrain (equation (18) in Luscher2017)
  // Lagrangian strain E_thermal = 1/2 (F_thermal^T F_thermal - I)
  // finite strain formula (Lubarda2002): F_thermal = exp((alpha/3)*(T-T_ref)) I
  thermal_eigenstrain = (1.0 / 2.0)
                      * (std::exp((2.0/3.0) * _thermal_expansion * (temp - _reference_temperature)) - 1.0)
                      * iden;

  // deviatoric stress + damage (equation (18) in Luscher2017): C : (Ee - alpha)
  pk2_new += xfac * _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);

  // Undamaged second piola-kirchoff stress to calculate undamaged plastic work
  _pk2_undamaged[_qp] += _elasticity_tensor[_qp] * (ee - thermal_eigenstrain);

  Real delta;
  // Pcor = correcting pressure = linearized form of the EOS
  // equation (18) in Luscher2017
  delta = 1.5 * (std::pow(Je , 2.0/3.0) - 1.0);
  pk2_new -= xfac * Kb * std::pow(Je , 2.0/3.0)
           * (delta * iden - 3.0 * thermal_eigenstrain)
           * invce;

  // Undamaged second piola-kirchoff stress to calculate undamaged plastic work
  _pk2_undamaged[_qp] -= Kb * std::pow(Je , 2.0/3.0)
                       * (delta * iden - 3.0 * thermal_eigenstrain)
                       * invce;

  // volumetric free energy = Psi_EOS in Luscher2017
  if (Je > 1.0) { // only in expansion
    _W0e[_qp] = (9.0 / 8.0) * Kb - (9.0 / 16.0) * Kb * (KT0prime - 4.0)
              + (9.0 / 8.0) * Kb * (1.0 - (3.0 / 2.0) * (KT0prime - 4.0)) * std::pow(V0V , 4.0/3.0)
              - (9.0 / 4.0) * Kb * (1.0 - (3.0 / 4.0) * (KT0prime - 4.0)) * std::pow(V0V , 2.0/3.0)
              + (9.0 / 16.0) * Kb * (KT0prime - 4.0) * std::pow(V0V , 2.0)
              + (Kb / _n_Murnaghan)
              * (std::exp(-_n_Murnaghan * _thermal_expansion * (temp - _reference_temperature)) - 1.0)
              * (Je - 1.0);
  } else {
    _W0e[_qp] = 0.0;
  }

  // volumetric coupling free energy = Psi_cpl in Luscher2017
  RankTwoTensor elastic_energy_tensor, thermal_coupling_tensor;

  elastic_energy_tensor = _elasticity_tensor[_qp] * ee;
  elastic_energy_tensor = 0.5 * ee * elastic_energy_tensor;
  _W0e[_qp] += elastic_energy_tensor.trace(); // 1/2 * Ee : C : Ee

  thermal_coupling_tensor = _elasticity_tensor[_qp] * thermal_eigenstrain; // C : alpha in equation 15 of Luscher2017
  thermal_coupling_tensor = ee * thermal_coupling_tensor;
  _W0e[_qp] -= thermal_coupling_tensor.trace(); // - Ee : C : alpha in equation 15 of Luscher2017

  _W0e[_qp] -= 0.5 * Kb * delta * delta; // - 1/2 K delta^2 in equation 15 of Luscher2017

  _W0e[_qp] += Kb * delta * thermal_eigenstrain.trace(); // + K delta alpha_v in equation 15 of Luscher2017

  // Used in PFFracBulkRate Jacobian
  // approximation: it should be only the positive part, undamaged pk2 is all
  _dW0e_dstrain[_qp] = _pk2_undamaged[_qp];

  // Used in StressDivergencePFFracTensors Jacobian
  // same approximation as above
  if (c < 1.0)
    _dstress_dc[_qp] = -_pk2_undamaged[_qp] * (2.0 * (1.0 - c));
  else
    _dstress_dc[_qp].zero();

  // Calculate bulk viscosity damping
  // C0 * dot(J) / J * |dot(J) / J| + C1 * dot(J) / J
  // C0 should be chosen of the order of rho * Le^2, rho = density, Le = element size
  // C1 should be chosen of the order of rho * Le * cs, cs = sound speed
  // Maheo et al. Mechanics Research Communications 38 (2011) 81 88
  trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  trD /= _deformation_gradient_old[_qp].det();

  pk2_new.addIa( _C0 * trD * std::abs(trD) );
  pk2_new.addIa( _C1 * trD );

  resid = _pk2_tmp - pk2_new;
}

// Calculate slip increment,dslipdtau. Override to modify.
void
FiniteStrainCrystalPlasticityDamageAmorEOS::getSlipIncrements()
{
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _slip_incr(i) = _a0(i) * std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i)) *
                    copysign(1.0, _tau(i)) * _dt;
    if (std::abs(_slip_incr(i)) > _slip_incr_tol * _dt)
    {
      _slip_incr(i) = _slip_incr_tol * _dt * copysign(1.0, _tau(i));
      //_err_tol = true;
#ifdef DEBUG
      mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
#endif
      return;
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dslipdtau(i) = _a0(i) / _xm(i) *
                    std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i) - 1.0) / _gss_tmp[i] *
                    _dt;
    if (std::abs(_slip_incr(i)) > _slip_incr_tol * _dt)
    {
      _dslipdtau(i) = 0.0;
    }
  }
}
