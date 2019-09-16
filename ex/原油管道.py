import numpy as np
"------------------------------------------Initialization--------------------------------------------------------------"
G = 4750 * np.e**7 / (350 * 24 * 3600) # 油品的质量流量
e = 5 * 10 ** (-5)
g = 9.8

'待确定参数'
t = 0
mu20 = 0
Q = 0
d = 1
nu = 1
L = 1
A = 1
l = 1
m = 1
Hsz = 1
Hc = 1
i = 1
hm = 1

epsilon_t = 1.825 - 0.001315 * mu20
mu0 = mu20 - epsilon_t * (t-2)
Q0 = G / mu0

Re = 4 * Q /(np.pi * d * nu)
epsilon = 2 * e / d
Re1 = 59.5 / (epsilon ** (8/7))

Re2 = (655 - 765 * np.log(epsilon))/epsilon

if Re <= 2000:
    lanbda = 64 / Re
if 2000 < Re and Re <= Re1:
    lanbda = 0.3164 * Re ** (-0.25)
if Re1 < Re and Re <= Re2:
    lanbda = 0.11 * (e/d + 68/Re) ** 0.25
if Re > Re2:
    lanbda = 1 / (1.74 - 2 * np.log(epsilon)) ** 2

'沿程摩阻损失'
upsilon = 4 * Q / (np.pi * d ** 2)
h_l = lanbda * L / d * upsilon ** 2 / (2 * g)
'综合摩阻损失'
beta = 8 * A / (4 ** m * np.pi ** (2-m) * g)
h_l = beta * l * Q ** (2-m) * upsilon ** m / d ** (5 - m)

H = h_l + Hsz

N= (i * L + Hsz) / Hc - hm

"------------------------------------------compute_temperature---------------------------------------------------------"
'待确定参数'
K = 1  # 管道总传热系数，W/(m2*°C)
D = 1  # 管道外径，m
c = 1  # 输油平均温度下油品的比热容，J/(kg*°C)
TR = 40  # 管道起点油温，°C
T0 = 20  # 周围介质温度，°C

'距离起始点l产生的温降'
'a,b参数'
a = K * np.pi * D / (G * c)
b = g * i / c / a

TL = (TR - T0 - b) / np.exp(a * L) + T0 + b