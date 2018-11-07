from dolfin import *
import matplotlib.pyplot as plt
import subprocess

set_log_level(30)
# Constants for the problem
epsilon_1 = Constant(3.5/(24*3600))
kappa = Constant(0.285)
M_s1 = Constant(3.6)
A_1hat = Constant(1)
a_T = Constant(0.4593)
Q = Constant(5)
# Define mesh
mesh = IntervalMesh(10, 0, 6)

# Define function space
P1 = FiniteElement('CG', interval, 5)
FS = FunctionSpace(mesh, P1)

# Define trial and test functions
u = Function(FS)
u_n = Function(FS)
f_1 = TestFunction(FS)

# Constants for the program
T = 10     # final time
num_steps = 500    # number of time steps
dt = T / num_steps
k = Constant(dt)

# Initial conditions
u_0 = Expression("cos(x[0])", degree = 1)
u_n = interpolate(u_0, FS)

# Define variational problem
F = ((u - u_n) / k)*f_1*dx - grad(u)[0]*f_1*dx
 
# Create VTK files for visualization output
vtkfile_v = File('reaction_system/velocity.pvd')
vtkfile_T = File('reaction_system/temperature.pvd')

# Solve the system for each time step
t = 0
for n in range(num_steps):
	if n%50 == 0:
		# Plot solution
		plot(u_n)
		plt.show()
		print("blaf", n)
		print(u_n(1))
	J = derivative(F, u)
	solve(F == 0, u, J=J)
	# Update previous solution
	u_n.assign(u)

p = subprocess.Popen('ls', stdout=subprocess.PIPE, shell=True)
print(p.communicate())