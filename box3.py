import dolfin as dol
import matplotlib.pyplot as plt
import plot
from const import *

# Set log-level to WARNING
dol.set_log_level(30)

# Time stepping
T = 5  # final time
num_steps = 200   # number of time steps
dt = T / num_steps
k = Constant(dt)

# Define mesh
mesh = dol.RectangleMesh(dol.Point(0,0), dol.Point(aspect_India, 1), 10, 10, "right/left")
"""
dol.plot(mesh)
plt.show()
#"""
# Define function space
P1 = dol.FiniteElement('CG', dol.triangle, 5)
element = dol.MixedElement([P1, P1, P1])
FS = dol.FunctionSpace(mesh, element)

# Define boundary classes
class Left(dol.SubDomain):
    def inside(self, x, on_boundary):
        return dol.near(x[0], 0.0)

class Right(dol.SubDomain):
    def inside(self, x, on_boundary):
        return dol.near(x[0], aspect_India)

class Bottom(dol.SubDomain):
    def inside(self, x, on_boundary):
        return dol.near(x[1], 0.0)

class Top(dol.SubDomain):
    def inside(self, x, on_boundary):
        return dol.near(x[1], 1.0)

# Initialize boundary instances
left = Left()
top = Top()
right = Right()
bottom = Bottom()

# Initialize mesh function for boundary domains
boundaries = dol.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)

# Boundary conditions
bcs = [dol.DirichletBC(FS, dol.Expression(('0','0', '0'), degree = 2), boundaries, 2),
       dol.DirichletBC(FS, dol.Expression(('5','10', '1'), degree = 2), boundaries, 4),
       dol.DirichletBC(FS, dol.Expression(("5*(1-x[1])", "10*(1-x[1])","1-x[1]"), degree = 2), boundaries, 1),
       dol.DirichletBC(FS, dol.Expression(("5*(1-x[1])", "10*(1-x[1])","1-x[1]"), degree = 2), boundaries, 3)]


# Define test functions
f_1, f_2, f_3 = dol.TestFunctions(FS)

# Split system function to access components
U = dol.Function(FS)
u, v, T = dol.split(U)

# Initial conditions
U_n = dol.interpolate(dol.Expression(("5*(1-x[1])", "10*(1-x[1])","1-x[1]"), degree = 2), FS)
#
u_n, v_n, T_n = dol.split(U_n)
dx, grad, derivative, solve = dol.dx, dol.grad, dol.derivative, dol.solve
# Define variational problem
F = ((u - u_n) / k)*f_1*dx - f*v*f_1*dx + kappa*grad(T)[0]*f_1*dx + epsilon_1*u*f_1*dx +\
    ((v - v_n) / k)*f_2*dx + f*u*f_2*dx + kappa*grad(T)[1]*f_2*dx + epsilon_1*v*f_2*dx +\
    A_1hat*((T - T_n) / k)*f_3*dx + M_s1*(grad(u)[0]+grad(v)[1])*f_3*dx - Q*f_3*dx + a_T*(v*grad(T)[1]+u*grad(T)[0])*f_3*dx


# Solve the system for each time step
t = 0
pltr_T = plot.Plotter2(mesh, id_ = '')
pltr_vel = plot.Plotter2(mesh, id_ = '')
temp = lambda x,y: T_n([x,y])
vel = lambda x,y: (u_n([x,y]), v_n([x,y]))
for n in range(num_steps):
    t += dt
    if n%1 == 0:
        pltr_vel.stream(vel,'qtcm1/velocity/', n, t, quantity ='box3_vel')
        pltr_T.plot(temp,'qtcm1/velocity/', n, t, quantity ='box3_temp')
    # Solve variational problem for time step
    J = derivative(F, U)
    solve(F == 0, U, bcs, J = J)
    U_n.assign(U)

pltr_vel.create_video()
pltr_T.create_video()
#"""