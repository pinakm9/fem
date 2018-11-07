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
T = 20	# final time
num_steps = 200		# number of time steps
dt = T / num_steps
k = Constant(dt)

# Define mesh
mesh = IntervalMesh(50, 0, 6)

# Define function space
P1 = FiniteElement('CG', interval, 3)
element = MixedElement([P1, P1])
FS = FunctionSpace(mesh, element)

# Define test functions
f_1, f_2 = TestFunctions(FS)

# Split system function to access components
u = Function(FS)
v, T = split(u)

# Initial conditions
u_n = interpolate(Expression(("x[0]", "cos(x[0])"), degree = 1), FS)
v_n, T_n = split(u_n)

# Define variational problem
F_1 = ((v - v_n) / k)*f_1*dx + epsilon_1*v*f_1*dx + kappa*grad(T)[0]*f_1*dx +\
	A_1hat*((T - T_n) / k)*f_2*dx + a_T*v*grad(T)[0]*f_2*dx + M_s1*grad(v)[0]*f_2*dx - Q*f_2*dx

F_2 = ((v - v_n) / k)*f_1*dx + epsilon_1*v*f_1*dx + beta*v*grad(v)[0]*f_1*dx + kappa*grad(T)[0]*f_1*dx +\
	A_1hat*((T - T_n) / k)*f_2*dx + a_T*v*grad(T)[0]*f_2*dx + M_s1*grad(v)[0]*f_2*dx - Q*f_2*dx


# Create VTK files for visualization output
"""vtkfile_v = File('qtcm1/velocity.pvd')
vtkfile_T = File('qtcm1/temperature.pvd')"""
pltr = Plotter(mesh, 0, 6, 'y', 'temperature', 8)
# Solve the system for each time step
t = 0
for n in range(num_steps):
	pltr.plot(lambda y: T_n([y]),'qtcm1/temperature/', n, t)
	t += dt
	# Solve variational problem for time step
	J = derivative(F_2, u)
	solve(F_2 == 0, u, J = J)
	# Save solution to file (VTK)
	"""_v, _T = u.split()
	vtkfile_v << (_v, t)
	vtkfile_T << (_T, t)"""
	# Update previous solution
	u_n.assign(u)

pltr.create_video()