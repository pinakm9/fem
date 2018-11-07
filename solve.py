from dolfin import *
import matplotlib.pyplot as plt
from plot import *
# Set log-level to WARNING
set_log_level(30)

# Constants for the problem
epsilon_1 = Constant(3.5/(24*3600))
kappa = Constant(0.285)
M_s1 = Constant(3.6)
A_1hat = Constant(0.4593)
a_T = Constant(0.0314)
Q = Constant(0)
beta = Constant(0.1)
# Constants for the program
T = 1  # final time
num_steps = 50    # number of time steps
dt = T / num_steps
k = Constant(dt)
f = 0.01

# Define mesh
mesh = Mesh("mesh/propeller_2d_coarse.xml.gz")
#plot(mesh)
#plt.show()
# Define function space
P1 = FiniteElement('P', triangle, 3)
element = MixedElement([P1, P1, P1])
FS = FunctionSpace(mesh, element)

# Define test functions
f_1, f_2, f_3 = TestFunctions(FS)

# Split system function to access components
U = Function(FS)
u, v, T = split(U)

# Initial conditions
U_n = interpolate(Expression(("pow(x[0],2)+x[1]", "cos(x[0])+cos(x[1])", "sin(x[0]+x[1])"), degree = 1), FS)
u_n, v_n, T_n = split(U_n)

# Define variational problem
F = ((u - u_n) / k)*f_1*dx - f*v*f_1*dx + kappa*grad(T)[0]*f_1*dx - epsilon_1*u*f_1*dx +\
      ((v - v_n) / k)*f_2*dx - f*u*f_2*dx + kappa*grad(T)[1]*f_2*dx - epsilon_1*v*f_2*dx +\
       A_1hat*((T - T_n) / k)*f_3*dx + M_s1*(grad(u)[0]+grad(v)[1])*f_3*dx - Q*f_2*dx


# Solve the system for each time step
t = 0
pltr = Plotter2(mesh)
temp = lambda x,y: v_n([x,y])
for n in range(num_steps):
    t += dt
    if n%1 == 0:
        pltr.plot(temp,'qtcm1/temperature/', n, t, quantity ='u-velocity')
    # Solve variational problem for time step
    J = derivative(F, U)
    solve(F == 0, U, J = J)
    U_n.assign(U)

#pltr.create_video()