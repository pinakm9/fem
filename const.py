from dolfin import Constant

# Constants for the problem
epsilon_1 = Constant(3.5/(24*3600))
kappa = Constant(0.285)
M_s1 = Constant(3.6)
A_1hat = Constant(0.4593)
a_T = Constant(0.0314)
Q = Constant(5)
K_H = Constant(1)
f = -0.00014
aspect_India =  2933/3214.0