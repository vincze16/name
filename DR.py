# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 10:44:00 2022

@author: vincz
"""
#DIFERENCIÁLNE ROVNICE

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid


tmin = 0 #t_0 zo zadanie initial - value problému
tmax = 5 #koncoví bod intervalu, na ktorom hľadáme riešenie
N = 100
n_iter = 5
disp_iter = 5
disp_points = 10

t, dt = np.linspace(tmin, tmax, N + 1, retstep = True)

ys = np.empty((n_iter + 1, len(t)))
y0 = 1
def func(t, y):
    #return 1/(2*y)
    return (2*y)/(t + 1)
def hladana_funk(t):
    #return np.sqrt(t + 1)
    return (t + 1)**2


y = np.full_like(t, y0)
ys[0, :] = y.copy()
for i in range(n_iter):
    y = y0 + cumulative_trapezoid(func(t,y), dx = dt, initial = 0)
    ys[i + 1, :] = y.copy()

plt.figure()
plt.plot(t, ys[:disp_iter + 1, :].T)
plt.scatter(t[::disp_points], hladana_funk(t[::disp_points]), c = "k", marker = "X")
plt.grid()
plt.title("Niekoľko prvých Pickardových iterácií", fontsize = 25)
plt.tight_layout()
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
ax = plt.gca()
ax.set_xlim([tmin, tmax])

y_A = np.zeros((1, len(t)))[0]
for i in range(1, len(t)):
    y_A_ = ys[:,i]
    while len(y_A_) >= 3:
        y_A_ = aitken(y_A_)
    y_A[i] = y_A_[-1]
y_A[0] = ys[-1,:][0]

plt.scatter(t[::disp_points], y_A[::disp_points], c = 'r', marker = 'X')

#Odchýlky
error = np.zeros((3,disp_points + 1))
error[0] = t[::disp_points]
for i in range(disp_points + 1):
    error[1][i] = abs(hladana_funk(error[0][i]) - [ys[-1][p] for p in range (len(t)) if t[p] == error[0][i]][0])
    error[2][i] = abs(hladana_funk(error[0][i]) - [y_A[p] for p in range (len(t)) if t[p] == error[0][i]][0])
