import numpy as np
"------------------------------------------Initialization--------------------------------------------------------------"
from sympy import *
from sympy.abc import a, b

T0 = 40
h1 = 2
h2 = 3
λ = 0.269
x1 = 0
x2 = 14.4
t_1 = 40
t_2 = 20

aa = solve([h1 * (a * x1 + b - t_1) + λ * a, h2 * (a * x2 + b - t_2) + λ * a], [a, b])
c1 = aa[a]
c2 = aa[b]
print(c1)
print(c2)