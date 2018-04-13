/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef COMPUTEPOLYCRYSTALELASTICITYTENSORCP_H
#define COMPUTEPOLYCRYSTALELASTICITYTENSORCP_H

#include "ComputeElasticityTensorBase.h"
#include "GrainDataTracker.h"
#include "RankTwoTensor.h"
#include "EulerAngleProvider.h"


//Forward Declarations
class ComputePolycrystalElasticityTensorCP;
class EulerAngleProvider;

template<>
InputParameters validParams<ComputePolycrystalElasticityTensorCP>();

/**
 * Compute an evolving elasticity tensor coupled to a grain growth phase field model.
 */
class ComputePolycrystalElasticityTensorCP : public ComputeElasticityTensorBase
{
public:
  ComputePolycrystalElasticityTensorCP(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();
  Real _length_scale;
  Real _pressure_scale;
//  EulerAngleProvider & _euler;

  /// Grain tracker object
  const GrainDataTracker<RankFourTensor> & _grain_tracker;
  const GrainDataTracker<EulerAngles> & _grain_tracker_crysrot;

  /// Number of order parameters
  unsigned int _op_num;

  /// Order parameters
  std::vector<const VariableValue *> _vals;

  /// vector of elasticity tensor material properties
  std::vector< MaterialProperty<RankFourTensor> *> _D_elastic_tensor;
  MaterialProperty<RankTwoTensor> & _crysrot;
 // MaterialProperty<RealVectorValue> & _angle2;
    
//  const ElementPropertyReadFile * _read_prop_user_object;
    
 // MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;

 // RotationTensor _R;
  /// Conversion factor from J to eV
  const Real _JtoeV;
};

#endif //COMPUTEPOLYCRYSTALELASTICITYTENSORCP_H
