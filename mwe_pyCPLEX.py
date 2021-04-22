#!/usr/bin/env python3

import numpy as np
from pymiqcqp._cplexmiqcqp_wrapper import solveMIQCQP

# Objective
Q0 = np.array([[33.0, -6.0, 0.0], [-6.0, 22.0, -23.0/2], [0.0, -23.0/2, 11.0]])
c0 = np.array([-1.0, -2.0, -3.0])
r0 = 0.0

# Linear constraitns
A1 = np.array([[-1.0, 1.0, 1.0], [1.0, -3.0, 1.0]])
b1 = np.array([20.0, 30.0])
A2 = np.array([[0.0, 0.0, 0.0]])
b2 = np.array([0.0])

# Quadratic constraints
Qs = np.eye(3).reshape(1, 3, 3)
cs = np.zeros(3).reshape(1, 3)
rs = np.array([-1.0])

# Variable bounds
lbs = -1000000.0*np.ones(3)
ubs = 1000000.0*np.ones(3)

# Integer variables x0 and x1
integer = np.array([0, 1])


print(solveMIQCQP(Q0, c0, r0, Qs, cs, rs, A1, b1, A2, b2, integer, lbs, ubs))
