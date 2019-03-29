# 2d_block_compress_neohookean_abaqus.cae
using abaqus build in Neo-Hookean model solve plane strain 2D block compress problem. the CPE6MH element is used.

the displacement load is 0.5, at t=0.91 (strain almost at 46%) the system generate negative eigenvalue. 

# 2d_block_compress_linear_umat.cae
using UMAT(Isotropic_Isothermal_Elasticity_UMAT.for) solve plane strain 2D block compress problem.

result can convergence.

# 2d_block_compress_abaqus_build_in_neo_hooken_umat_compressible.cae
using UMAT(neo_hookean_abaqus_compressible_umat.for) solve plane strain 2D block compress problem.

result can convergence.

# 2d_block_compress_my_neo_hooken_umat_compressible.cae
using UMAT(neo_hookean_my_compressible_umat.for) solve plane strain 2D block compress problem.

result can convergence.

the displacement load is 0.5,at t=0.914 (strain almost 46%) the system generate negative eigenvalue.

the result is same with abaqus build in model