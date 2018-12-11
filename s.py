from dolfin import *
from plot import Plotter
import numpy as np


# Set log-level to WARNING
set_log_level(30)

# Constants for the problem
epsilon_1 = Constant(3.5/(24*3600))
kappa = Constant(0.285)
M_s1 = Constant(3.6)
A_1hat = Constant(0.4593)
a_T = Constant(0.0314)
Q = Constant(5)
beta = Constant(0.1)
# Constants for the program
T = 2	# final time
num_steps = 200		# number of time steps
dt = T / num_steps
k = Constant(dt)

# Define mesh
mesh = IntervalMesh(10, 0, 1)

# Define function space
P1 = FiniteElement('P', interval, 2)
element = MixedElement([P1, P1])
FS = FunctionSpace(mesh, element)
R = FunctionSpace(mesh, 'Real', 0)

# Define test functions

f_1, f_2 = TestFunctions(FS)

s = TestFunctions(FunctionSpace(mesh, 'Real', 0))
r = TrialFunctions(R)
# Split system function to access components
u = Function(FS)
v, T = split(u)

left = Expression("1-x[0]", degree = 1)
right = Expression("x[0]", degree = 1)

# Initial conditions
u_n = interpolate(Expression(("sin(3.14159265359*x[0])", "30"), degree = 2), FS)
v_n, T_n = split(u_n)

u_D = Expression(('0', '30'), degree = 2, t=0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(FS, u_D, boundary)


# Define variational problem
F_1 = ((v - v_n) / k)*f_1*dx + epsilon_1*v*f_1*dx + kappa*grad(T)[0]*f_1*dx +\
	A_1hat*((T - T_n) / k)*f_2*dx + a_T*grad(T)[0]*v*f_2*dx + M_s1*grad(v)[0]*f_2*dx - Q*f_2*dx #+ kappa*10*(right*f_1-left*f_1)*ds
 

F_2 = ((v - v_n) / k)*f_1*dx + epsilon_1*v*f_1*dx + beta*v*grad(v)[0]*f_1*dx + kappa*grad(T)[0]*f_1*dx +\
	A_1hat*((T - T_n) / k)*f_2*dx + a_T*v*grad(T)[0]*f_2*dx + M_s1*grad(v)[0]*f_2*dx - Q*f_2*dx

# Create VTK files for visualization output
"""vtkfile_v = File('qtcm1/velocity.pvd')
vtkfile_T = File('qtcm1/temperature.pvd')"""
pltr = Plotter(mesh, id_='4')
# Solve the system for each time step
t = 0
v_ = lambda y: v_n([y])
T_ = lambda y: T_n([y])
for n in range(num_steps):
	#pltr.plot(v_,'qtcm1/velocity/', n, t, quantity = 'velocity_42')
	pltr.plot(T_,'qtcm1/velocity/', n, t, quantity = 'temp_43')
	t += dt
	# Solve variational problem for time step
	J = derivative(F_1, u)
	solve(F_1 == 0, u,bc,  J = J)
	# Save solution to file (VTK)
	"""_v, _T = u.split()
	vtkfile_v << (_v, t)
	vtkfile_T << (_T, t)"""
	# Update previous solution
	u_n.assign(u)

pltr.create_video()