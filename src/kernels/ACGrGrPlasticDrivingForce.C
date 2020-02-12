//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACGrGrPlasticDrivingForce.h"

#include "Material.h"
#include "RankFourTensor.h"
#include "RankTwoTensor.h"

registerMooseObject("veryhappyApp", ACGrGrPlasticDrivingForce);

template <>
InputParameters
validParams<ACGrGrPlasticDrivingForce>()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Adds plastic energy and elastic energy contribution to the Allen-Cahn equation");
  params.addRequiredParam<MaterialPropertyName>(
      "D_tensor_name", "The elastic tensor derivative for the specific order parameter");
  params.addRequiredParam<MaterialPropertyName>(
      "D_energy_name", "The plastic energy derivative for the specific order parameter");
  return params;
}

ACGrGrPlasticDrivingForce::ACGrGrPlasticDrivingForce(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _D_plastic_energy(getMaterialProperty<Real>("D_energy_name")),
    _D_elastic_tensor(getMaterialProperty<RankFourTensor>("D_tensor_name")),
    _elastic_strain(getMaterialPropertyByName<RankTwoTensor>("elastic_strain"))
{
}

Real
ACGrGrPlasticDrivingForce::computeDFDOP(PFFunctionType type)
{
  // Access the heterogeneous strain calculated by the Solid Mechanics kernels
  RankTwoTensor strain(_elastic_strain[_qp]);

  // Compute the partial derivative of the stress wrt the order parameter
  RankTwoTensor D_stress = _D_elastic_tensor[_qp] * strain;
    
  Real D_pla_energy_tmp = 0.0;
   // D_pla_energy_tmp = 1.0;
  D_pla_energy_tmp = _D_plastic_energy[_qp];
    
  //mooseWarning("Kernel D_plastic_energy = " , D_pla_energy_tmp);
  //mooseWarning("Kernel D_elastic_energy = " , D_stress.doubleContraction(strain));
    
  switch (type)
  {
    case Residual:
      return 0.5 *
             D_stress.doubleContraction(strain) + D_pla_energy_tmp; // Compute the deformation energy driving force

    case Jacobian:
      return 0.0;
  }

  mooseError("Invalid type passed in");
}
