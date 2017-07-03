[Mesh]
  type = FileMesh
  file = hmxsylgard100500.msh
  second_order = true
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
    order = SECOND
    family = LAGRANGE
  [../]
  [./disp_y]
    order = SECOND
    family = LAGRANGE
  [../]
  [./b]
    order = SECOND
    family = LAGRANGE
  [../]
  [./c]
    order = SECOND
    family = LAGRANGE
  [../]
  [./temp]
    order = SECOND
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vel_x]
    order = SECOND
    family = LAGRANGE
  [../]
  [./accel_x]
    order = SECOND
    family = LAGRANGE
  [../]
  [./vel_y]
    order = SECOND
    family = LAGRANGE
  [../]
  [./accel_y]
    order = SECOND
    family = LAGRANGE
  [../]
[]

[ICs]
  [./cIC]
    type = FunctionIC
    variable = c
    block = 'hmx'
    function = xyfunc_c
  [../]
  [./bIC]
    type = FunctionIC
    variable = b
    block = 'hmx'
    function = xyfunc_b
  [../]
  [./v_ic_bnd]
    type = ConstantIC
    variable = vel_x
    boundary = left
    value = 0.1
  [../]
  [./tempIC]
    type = FunctionIC
    variable = temp
    block = 'hmx'
    function = tempfunc
  [../]
[]

[Functions]
  [./tfunc]
    type = ParsedFunction
    value = 0.1
  [../]
  [./xyfunc_c]
    type = ParsedFunction
    value = 'min(exp(-abs(y-250)/10)*if((x-150)/100,0,1)+exp(-abs(y-350)/10)*if((x-550)/100,0,1) +exp(-(sqrt(2)/10)*abs((x-y)/2-50))*max(0.0,floor(4.925-x/100))*(1.0-max(0.0,min(1.0,floor(3.925-x/100)))),1)'
  [../]
  [./xyfunc_b]
    type = ParsedFunction
    value = '(1/(10*10))*min(exp(-abs(y-250)/10)*if((x-150)/100,0,1)+exp(-abs(y-350)/10)*if((x-550)/100,0,1) +exp(-(sqrt(2)/10)*abs((x-y)/2-50))*max(0.0,floor(4.925-x/100))*(1.0-max(0.0,min(1.0,floor(3.925-x/100)))),1)'
  [../]
  [./tempfunc]
    type = ParsedFunction
    value = 293.0
  [../]
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = 0.3025
    gamma = 0.6
  [../]
  [./inertia_y]
    type = InertialForce
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = 0.3025
    gamma = 0.6
  [../]
  [./pfbulk]
    type = PFFracBulkRate
    variable = c
    l = 10
    beta = b
    visco = 100
    gc_prop_var = 'gc_prop'
    G0_var = 'G0_pos'
    dG0_dstrain_var = 'dG0_pos_dstrain'
  [../]
  [./DynamicTensorMechanics]
    displacements = 'disp_x disp_y'
    #use_displaced_mesh = true
  [../]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
  [../]
  [./dcdt]
    type = TimeDerivative
    variable = c
  [../]
  [./pfintvar]
    type = Reaction
    variable = b
  [../]
  [./pfintcoupled]
    type = PFFracCoupledInterface
    variable = b
    c = c
  [../]
  [./hc]
    type = HeatConduction
    variable = temp
    diffusion_coefficient = thermal_conductivity
  [../]
  [./hct]
    type = HeatConductionTimeDerivative
    variable = temp
    specific_heat = specific_heat
    density_name = density
  [../]
  [./thermalexpansionheat]
    type = CrackPropagationHeatSource
    variable = temp
    c = c
    G0_var = 'G0_pos'
    dG0_dstrain_var = 'dG0_pos_dstrain'
    #block = 'hmx'
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
  [../]
  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
    index_j = 0
    index_i = 1
  [../]
  [./strain_xx]
    type = RankTwoAux
    variable = strain_xx
    rank_two_tensor = total_strain
    index_j = 0
    index_i = 0
  [../]
  [./strain_yy]
    type = RankTwoAux
    variable = strain_yy
    rank_two_tensor = total_strain
    index_j = 1
    index_i = 1
  [../]
  [./strain_xy]
    type = RankTwoAux
    variable = strain_xy
    rank_two_tensor = total_strain
    index_j = 0
    index_i = 1
  [../]
  [./accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = 0.3025
  [../]
  [./vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = 0.6
  [../]
  [./accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = 0.3025
  [../]
  [./vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = 0.6
  [../]
[]

[BCs]
  [./xdispleft]
    type = PresetVelocity
    variable = disp_x
    boundary = left
    function = tfunc
    velocity = 1.0
  [../]
  [./xdispright]
    type = PresetBC
    variable = disp_x
    boundary = right
    value = 0
  [../]
  [./ydispbottom]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./tright]
    type = DirichletBC
    variable = temp
    boundary = right
    value = 293
  [../]
[]

[Materials]
  [./elastic]
    type = LinearIsoElasticPFDamage
    c = c
    kdamage = 1e-6
  [../]
  [./pfbulkmatsylgard]
    type = PFFracBulkRateMaterial
    block = 'sylgard'
    gc = 0.4
  [../]
  [./pfbulkmathmx]
    type = PFFracBulkRateMaterial
    block = 'hmx'
    gc = 0.002
  [../]
  [./elasticity_tensor_sylgard]
    youngs_modulus = 0.159 #GPa
    poissons_ratio = 0.48
    type = ComputeIsotropicElasticityTensor
    block = 'sylgard'
  [../]
  [./elasticity_tensor_hmx]
    type = ComputeElasticityTensorCP
    block = 'hmx'
    C_ijkl = '22.2 9.6 13.2 0.0 0.1 0.0 23.9 13.0 0.0 -4.7 0.0 23.4 0.0 -1.6 0.0 9.2 0.0 -2.5 11.1 0.0 10.1'
    fill_method = symmetric21
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 'sylgard hmx'
    displacements = 'disp_x disp_y'
  [../]
  [./density_sylgard]
    type = GenericConstantMaterial
    block = 'sylgard'
    prop_names = 'density'
    prop_values = '1.03'
  [../]
  [./density_hmx]
    type = GenericConstantMaterial
    block = 'hmx'
    prop_names = 'density'
    prop_values = '1.91'
  [../]
  [./thermal_conductivity_hmx]
    type = GenericConstantMaterial
    block = 'hmx'
    prop_names = 'thermal_conductivity'
    prop_values = '0.8e-6'
  [../]
  [./thermal_conductivity_sylgard]
    type = GenericConstantMaterial
    block = 'sylgard'
    prop_names = 'thermal_conductivity'
    prop_values = '0.27e-6'
  [../]
  [./specific_heat_hmx]
    type = GenericConstantMaterial
    block = 'hmx'
    prop_names = 'specific_heat'
    prop_values = '995e-6'
  [../]
  [./specific_heat_sylgard]
    type = GenericConstantMaterial
    block = 'sylgard'
    prop_names = 'specific_heat'
    prop_values = '1465e-6'
  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201 hypre boomeramg 20'
  line_search = 'none'

#  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
#  petsc_options_value = 'lu    superlu_dist'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  #l_max_its = 20
  #nl_max_its = 500

  dt = 0.05
  dtmin = 1e-9
  num_steps = 200000
[]

[Outputs]
  exodus = true
  interval = 400
[]
