# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 19:40:38 2022

@author: vincz
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin

a = 1
b = 100
def rosenbrock(y):
    return (a - y[0])**2 + b*(y[1] - y[0]**2)**2

x0 = [-1, 2]
minimum = [a, a**2]

area = [0]*2
area[0] = [-2, 2]
area[1] = [-2, 3]
delenie = 501
oblast = np.meshgrid(np.linspace(area[0][0], area[0][1], delenie),np.linspace(area[1][0], area[1][1], delenie))
values = np.meshgrid(np.linspace(area[0][0], area[0][1], delenie),np.linspace(area[1][0], area[1][1], delenie))[0]
for i in range(delenie):
    for p in range(delenie):
        x_ = np.array([oblast[0][i][p], oblast[1][i][p]])
        values[i][p] = rosenbrock(x_)


N = 500 #10
chyba = [0]*4


#Cauchy
def lambda_f(t):
    return rosenbrock(x + t*g_k)

x = x0
postupnost = [x]

plt.figure()
plt.scatter(postupnost[-1][0], postupnost[-1][1], c = "r", marker = "X")
plt.contour(oblast[0], oblast[1], values, levels = 100)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)

for i in range(N):
    g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    Lambda_k = fmin(lambda_f, 0)
    x = x + Lambda_k*g_k
    postupnost.append(x)
    
    plt.plot([postupnost[-2][0], postupnost[-1][0]],[postupnost[-2][1], postupnost[-1][1]], c = "b")
    plt.scatter(postupnost[-1][0], postupnost[-1][1], c = "r", marker = "X")
plt.scatter(minimum[0], minimum[1], c = "m", marker = "X")
plt.gca().set_aspect('equal', adjustable = 'box')
chyba[0] = np.linalg.norm(postupnost[-1] - minimum)


#Konštantná dĺžka kroku
Lambda = 0.002 #0.002

x = x0
postupnost2 = [x]

plt.figure()
plt.scatter(minimum[0], minimum[1], c = "m", marker = "X")
plt.scatter(postupnost2[-1][0], postupnost2[-1][1], c = "r", marker = "X")
plt.contour(oblast[0],oblast[1], values, levels = 100)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)

for i in range(N):
    g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    x = x - Lambda*g_k
    postupnost2.append(x)
    
    plt.plot([postupnost2[-2][0], postupnost2[-1][0]],[postupnost2[-2][1], postupnost2[-1][1]], c = "b")
    plt.scatter(postupnost2[-1][0], postupnost2[-1][1], c = "r", marker = "X")
plt.scatter(minimum[0], minimum[1], c = "m", marker = "X")
plt.gca().set_aspect('equal', adjustable = 'box')
chyba[1] = np.linalg.norm(postupnost2[-1] - minimum)


#Heavy ball method
beta = 0.8

x = x0
postupnost3 = [x]

plt.figure()
plt.scatter(minimum[0], minimum[1], c = "m", marker = "X")
plt.scatter(postupnost3[-1][0], postupnost3[-1][1], c = "r", marker = "X")
plt.contour(oblast[0],oblast[1], values, levels = 100)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
x = x - Lambda*g_k
postupnost3.append(x)
plt.scatter(postupnost3[-1][0], postupnost3[-1][1], c = "r", marker = "X")
plt.plot([postupnost3[-2][0], postupnost3[-1][0]],[postupnost3[-2][1], postupnost3[-1][1]], c = "b")
for i in range(N-1):
    g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    x = x - Lambda*g_k + beta*(postupnost3[-1] - postupnost3[-2])
    postupnost3.append(x)
    plt.plot([postupnost3[-2][0], postupnost3[-1][0]],[postupnost3[-2][1], postupnost3[-1][1]], c = "b")
    plt.scatter(postupnost3[-1][0], postupnost3[-1][1], c = "r", marker = "X")
plt.scatter(minimum[0], minimum[1], c = "m", marker = "X")
plt.gca().set_aspect('equal', adjustable = 'box')
chyba[2] = np.linalg.norm(postupnost3[-1] - minimum)


#Barzilai&Borwein
x = x0

postupnost4 = [x]
plt.figure()
plt.scatter(minimum[0], minimum[1], c = "m", marker = "X")
plt.scatter(postupnost4[-1][0], postupnost4[-1][1], c = "r", marker = "X")
plt.contour(oblast[0],oblast[1], values, levels = 100)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)

g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
Lambda_k = fmin(lambda_f, 0)
x = x + Lambda_k*g_k
postupnost4.append(x)
plt.plot([postupnost4[-2][0], postupnost4[-1][0]],[postupnost4[-2][1], postupnost4[-1][1]], c = "b")

s_k = Lambda_k*g_k
g_k_ = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
y_k = g_k_ - g_k
alpha_k = np.dot(s_k.T, y_k)/np.dot(s_k.T,s_k)

for i in range(N-1):
    g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    s_k = -1/alpha_k*g_k
    x = x + s_k
    postupnost4.append(x)
    g_k_ = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    y_k = g_k_ - g_k
    alpha_k = np.dot(s_k.T, y_k)/np.dot(s_k.T,s_k)
    plt.plot([postupnost4[-2][0], postupnost4[-1][0]],[postupnost4[-2][1], postupnost4[-1][1]], c = "b")
    plt.scatter(postupnost4[-1][0], postupnost4[-1][1], c = "r", marker = "X")
plt.scatter(minimum[0], minimum[1], c = "m", marker = "X")
plt.gca().set_aspect('equal', adjustable = 'box')
chyba[3] = np.linalg.norm(postupnost4[-1] - minimum)



#####Potrebný počet iterácií#####

eps = 10**(-6)
potrebne_iter = [0]*4

#Cauchy
x = x0

postupnost_ = [x]
g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
while np.linalg.norm(g_k) > eps:
    g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    Lambda_k = fmin(lambda_f, 0)
    x = x + Lambda_k*g_k
    postupnost_.append(x)
potrebne_iter[0] = len(postupnost_) - 1


#Konštantná dĺžka kroku
x = x0
postupnost2_ = [x]
g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
while np.linalg.norm(g_k) > eps:
    g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    x = x - Lambda*g_k
    postupnost2_.append(x)
potrebne_iter[1] = len(postupnost2_) - 1


#Heavy ball
beta = 0.91

x = x0
postupnost3_ = [x]
g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
x = x - Lambda*g_k
postupnost3_.append(x)
g_k = np.dot(G,x) + h
while np.linalg.norm(g_k) > eps:
    g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    x = x - Lambda*g_k + beta*(postupnost3_[-1] - postupnost3_[-2])
    postupnost3_.append(x)
potrebne_iter[2] = len(postupnost3_) - 1


#Barzilai&Borwein
x = x0
postupnost4_ = [x]
g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
Lambda_k = fmin(lambda_f, 0)
x = x + Lambda_k*g_k
postupnost4.append(x)
plt.plot([postupnost4[-2][0], postupnost4[-1][0]],[postupnost4[-2][1], postupnost4[-1][1]], c = "b")
s_k = -Lambda_k*g_k
g_k_ = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
y_k = g_k_ - g_k
alpha_k = np.dot(s_k.T, y_k)/np.dot(s_k.T,s_k)
while np.linalg.norm(g_k) > eps:
    g_k = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    s_k = -1/alpha_k*g_k
    x = x + s_k
    postupnost4.append(x)
    g_k_ = np.array([-2*(a - x[0]) - 4*b*x[0]*(x[1] - x[0]**2), 2*b*(x[1] - x[0]**2)])
    y_k = g_k_ - g_k
    alpha_k = np.dot(s_k.T, y_k)/np.dot(s_k.T,s_k)
potrebne_iter[3] = len(postupnost4_) - 1