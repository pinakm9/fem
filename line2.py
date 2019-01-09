import dolfin as dol
import matplotlib.pyplot as plt
import plot
from const import *

# Set log-level to WARNING
dol.set_log_level(30)

# Time stepping
T = 3  # final time
num_steps = 300   # number of time steps
dt = T / num_steps
k = Constant(dt)

# Define mesh
mesh = dol.IntervalMesh(20, 0, 1)
"""
dol.plot(mesh)
plt.show()
#"""
# Define function space
P1 = dol.FiniteElement('CG', dol.interval, 2)
element = dol.MixedElement([P1, P1])
FS = dol.FunctionSpace(mesh, element)

# Define boundary classes
class Left(dol.SubDomain):
    def inside(self, x, on_boundary):
        return dol.near(x[0], 0.0)

class Right(dol.SubDomain):
    def inside(self, x, on_boundary):
        return dol.near(x[0], 1.0)

# Initialize boundary instances
left = Left()
right = Right()

# Initialize mesh function for boundary domains
boundaries = dol.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
right.mark(boundaries, 2)

# # Boundary conditions
bcs = [dol.DirichletBC(FS, dol.Expression(('5', '1'), degree = 2), boundaries, 1),
       dol.DirichletBC(FS, dol.Expression(('1', '2-exp(1)'), degree = 2), boundaries, 2)]
"""
bc_l = dol.DirichletBC(FS, dol.Expression(("5","0"), degree = 2), 'x[0] < DOLFIN_EPS')
bc_r = dol.DirichletBC(FS, dol.Expression(("0","1"), degree = 2), 'x[0] > 1- DOLFIN_EPS')
bc_left = dol.DirichletBC(FS, dol.Expression(("5","0"), degree = 2),left_boundary, method='pointwise')
bc_right = dol.DirichletBC(FS, dol.Expression(("0","1"), degree = 2),right_boundary, method='pointwise')

def boundary(x, on_boundary):
    return (abs(x[0]) < 1e-15) or (abs(x[0]-1.0)< 1e-15)
bc = dol.DirichletBC(FS, dol.Expression(("5*(1-x[0])","1-x[0]"), degree = 2), boundary)
"""
# Define test functions
f_1, f_2 = dol.TestFunctions(FS)

# Split system function to access components
U = dol.Function(FS)
v, T = dol.split(U)

# Initial conditions
U_n = dol.interpolate(dol.Expression(("5*(1-x[0])","2-exp(x[0])"), degree = 2), FS)
v_n, T_n = dol.split(U_n)
dx, grad, derivative, solve = dol.dx, dol.grad, dol.derivative, dol.solve
# Define variational problem
F = ((v - v_n) / k)*f_1*dx + epsilon_1*v*f_1*dx + kappa*grad(T)[0]*f_1*dx +\
    A_1hat*((T - T_n) / k)*f_2*dx + a_T*grad(T)[0]*v*f_2*dx + M_s1*grad(v)[0]*f_2*dx - Q*f_2*dx

L = kappa*f_1*dol.ds+M_s1*5*f_2*dol.ds
# Solve the system for each time step
t = 0
pltr_T = plot.Plotter(mesh, id_ = '')
pltr_vel = plot.Plotter(mesh, id_ = '')
temp = lambda y: T_n([y])
vel = lambda y: v_n([y])
for n in range(num_steps):
    t += dt
    if n%1 == 0:
        pltr_vel.plot(vel,'qtcm1/velocity/', n, t, quantity ='line2_vel')
        pltr_T.plot(temp,'qtcm1/velocity/', n, t, quantity ='line2_temp')
    # Solve variational problem for time step
    J = derivative(F, U)
    solve(F == 0, U, bcs, J = J)
    U_n.assign(U)

pltr_vel.create_video()
pltr_T.create_video()
#"""