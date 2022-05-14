# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 13:36:17 2022

@author: vincz
"""
#TAYLOROVE A FOURIEROVE RADY

from sympy import *
import matplotlib.pyplot as plt

#Taylorove rady
M = 6
f = E**x
margins = [-3, 3]

F = []
for i in range(M):
    F.append(series(f, x, x0 = 0, n = i + 1).removeO())


F_ = [F[i] for i in range (M-1) if F[i] != F[i + 1]]
if F_[-1] != F[-1]:
    F_.append(F[-1])
F = F_

while len(F_)>=3:
    F_ = aitken(F_)
for i in range(len(F_)):
    F_[i] = simplify(F_[i])
    
#Grafy
num_of_points = 1000
xlim = [-3, 3]
os_x = np.linspace(float(margins[0]), float(margins[1]), num_of_points)
f1 = np.zeros(num_of_points)
f2 = np.zeros(num_of_points)
f3 = np.zeros(num_of_points)
for i in range(num_of_points):
    f1[i] = f.evalf(20, subs = {x:os_x[i]})
    f2[i] = F[-1].evalf(20, subs = {x:os_x[i]})
    f3[i] = F_[-1].evalf(20, subs = {x:os_x[i]})
plt.figure()
plt.plot(os_x, f1)
plt.plot(os_x, f2)
plt.plot(os_x, f3)
plt.xlabel('x', fontsize = 25)
plt.ylabel('f(x)', fontsize = 25)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
ax = plt.gca()
ax.set_xlim(xlim)
plt.grid()

#Kontrola zlepšenia presnosti
test_bod = pi
des_miesta = 20
f_T = f - F[-1]
f_A = f - F_[-1]
chyba_T = [0]*2
chyba_T[0] = abs(f_T).evalf(des_miesta, subs = {x:test_bod})
chyba_T[1] = abs(f_A).evalf(des_miesta, subs = {x:test_bod})


#Fourierove rady
P = 5
g = x
b1 = -pi
b2 = pi
margins_ = [-2*pi, 2*pi]

a_ = [0]*(P + 1)
b_ = [0]*(P + 1)
S = [0]*(P + 1)
for i in range (P + 1):
    a_[i] = 2/(b2 - b1)*integrate(g*cos(2*pi/(b2 - b1)*i*x), (x, b1, b2))
    b_[i] = 2/(b2 - b1)*integrate(g*sin(2*pi/(b2 - b1)*i*x), (x, b1, b2))
S[0] = a_[0]/2
for i in range(1, P + 1):
    S[i] = S[i - 1] + a_[i]*cos(2*pi/(b2 - b1)*i*x) + b_[i]*sin(2*pi/(b2 - b1)*i*x)

S_ = [S[i] for i in range (P) if S[i] != S[i + 1]]
if S_[-1] != S[-1]:
    S_.append(S[-1])
S = S_

while len(S_) >= 3:
    S_ = aitken(S_)
S_[-1] = simplify(S_[-1])

#Kontrola zlepšenia presnosti
test_bod = -0.25
des_miesta = 20
f_F = g-S_[-1]
f_S = g-S__[-1]
chyba_F = [0]*2
chyba_F[0] = abs(f_T).evalf(des_miesta, subs={x:test_bod})
chyba_F[1] = abs(f_S).evalf(des_miesta, subs={x:test_bod})

#Grafy
num_of_points_ = 2000
xlim = [-6, 6]
ylim = [-6, 6]
os_x = np.linspace(float(margins_[0]), float(margins_[1]), num_of_points_)
f1_ = np.zeros(num_of_points_)
f2_ = np.zeros(num_of_points_)
f3_ = np.zeros(num_of_points_)
for i in range(num_of_points_):
    f1_[i] = g.evalf(20, subs = {x:os_x[i]})
    f2_[i] = S[-1].evalf(20, subs = {x:os_x[i]})
    f3_[i] = S_[-1].evalf(20, subs = {x:os_x[i]})
plt.figure()
plt.plot(os_x, f1_)
plt.plot(os_x, f2_)
plt.plot(os_x, f3_)
plt.xlabel('x', fontsize = 25)
plt.ylabel('f(x)', fontsize = 25)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
ax = plt.gca()
ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.grid()