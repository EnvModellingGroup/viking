# Simple rectangle with sigmoidal slope tsunami simulation
# ================================================================================
from thetis import *
from firedrake.petsc import PETSc
from math import sin,radians
import sys

def set_tsunami_field(elev, t):

    amplitude = -1.0
    mesh2d = elev.function_space().mesh()
    xvector = mesh2d.coordinates.dat.data
    evector = elev.dat.data
    for i,xy in enumerate(xvector):
        evector[i] = amplitude * sin(radians(t/15.))

    return evector

# the mesh. 

# 100 x 50 km mesh with 250x125 triangles
# scale this appropriately for the number of cores you want to use
# 250x125 is ~4-8 cores
mesh2d = RectangleMesh(250, 125, 100000., 50000.)

# we set the lefthand x edge as the incoming wave
physID = 1
#timestepping options
dt = 10 # reduce this if solver does not converge
t_export = 60.0
t_end = 7200. # 2 hours
output_directory = "output"
# model options
wd_alpha = 1.0
h_viscosity = 10.
manning_drag_coefficient = 0.025
use_wetting_and_drying = True

output_dir = create_directory(output_directory)

# Bathymetry and viscosity field
P1_2d = get_functionspace(mesh2d, 'DG', 1)
bathymetry_2d = Function(P1_2d, name='Bathymetry')
viscosity_2d = Function(P1_2d, name='viscosity')
# Bathymetry combination of bell curve and sloping bathymetry
x, y = SpatialCoordinate(mesh2d)
shift = 50000.
scale = 10000.
const = -510.
shift2 = 500.
# sigmoidal bathy function
bathymetry_2d.interpolate(const * (exp((x-shift)/scale) / (exp((x-shift)/scale) + 1. )) + shift2)

# viscosity, create a distance function
L = Constant(1e2)
V = FunctionSpace(mesh2d, 'CG', 1)
# Calculate distance to open boundary
bcs = DirichletBC(V, 0.0, physID)
v = TestFunction(V)
u = Function(V)#.interpolate(0.0)
u.interpolate(Constant(0.0))
solver_parameters={'snes_rtol': 1e-3,
         'ksp_type': 'preonly',
                   } 
# Before we solve the Eikonal equation, let's solve a Laplace equation to
# generate an initial guess
F = L**2*(inner(grad(u), grad(v))) * dx - v * dx

solve(F == 0, u, bcs, solver_parameters=solver_parameters)
chk = DumbCheckpoint('dist', mode=FILE_CREATE)
with timed_stage('initialising dist'):
    File('dist.pvd').write(u)
# epss values set the accuracy (in meters) of the final "distance to boundary" function. To make
# more accurate add in extra iterations, eg, 500., 250., etc. This may result in the solver not
# converging.
epss = [10000., 1000., 500., 250.]
# solve Eikonal equations
for i, eps in enumerate(epss):
    PETSc.Sys.Print('Solving Eikonal with eps == ' + str(float(eps)))
    F = inner(sqrt(inner(grad(u), grad(u))), v) * dx - v * dx + Constant(eps)*inner(grad(u), grad(v)) * dx
    solve(F == 0, u, bcs, solver_parameters=solver_parameters)

viscosity_2d.interpolate(max_value(h_viscosity, 1000 * (1. - u / 10e3)))

# --- create solver ---
solverObj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
options = solverObj.options
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = output_dir
options.fields_to_export = ['uv_2d','elev_2d']
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
# spatial discretisation
options.element_family = "dg-dg"
options.swe_timestepper_type='CrankNicolson'
options.manning_drag_coefficient = Constant(manning_drag_coefficient)
options.horizontal_viscosity = viscosity_2d
options.timestep = dt
options.use_wetting_and_drying = True
options.wetting_and_drying_alpha = Constant(wd_alpha)

# boundary conditions
tsunami_elev = Function(FunctionSpace(mesh2d, "CG", 1), name='tidal_elev')
solverObj.bnd_functions['shallow_water'] = {
    physID: {'elev': tsunami_elev}
    # any other boundary has u_n = 0.0
}

# Set up usual SWE terms
solverObj.create_equations()

def update_forcings(t):
    set_tsunami_field(tsunami_elev, t)        

update_forcings(0.0)
solverObj.assign_initial_conditions(uv=Constant((1e-9,0.0)), elev=Constant(1e-9))
# Run model
solverObj.iterate(update_forcings=update_forcings)


