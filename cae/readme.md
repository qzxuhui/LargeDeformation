2d_block_compress_neohookean_abaqus.cae
using abaqus build in Neo-Hookean model solve plane strain 2D block compress problem. the CPE6MH element is used.
the displacement load is 0.5, at t=0.91 (strain almost at 46%) the system generate negative eigen value. 

2d_block_compress_linear_umat.cae
using UMAT(Isotropic_Isothermal_Elasticity_UMAT.for) solve plane strain 2D block compress problem.

