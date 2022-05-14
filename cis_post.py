# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 19:37:55 2022

@author: vincz
"""
#ČÍSELNÉ POSTUPNOSTI

from sympy import *
import numpy as np
import matplotlib.pyplot as plt


#Definovanie postupnosti, ktorú sčítavame
def a_n(n):
    #return (-1)**n/(1 + 2*n)
    #return (-1)**n/2**n #Shanks nespolupracuje
    #return (-1)**n/(n + 1)
    #return (-1)**n/(n + 1)**2
    return n*(-2)**(-n)
    #return n/(2**n)
    #return (-2)**(-n)
    #return (n + 1)/(2**n)

#Definovanie počtu sčítavaných členov a spočítanie čiastočných súčtov
N = 8
a = [0]*N
for i in range(N):
    a[i] = a_n(i)
sums = [0]*N
sums[0] = a[0]
for i in range(1,N):
    sums[i] = sums[i - 1] + a[i]

#Skutočný súčet
SUM = Sum(a_n(x), (x, 0, oo)).doit()


#POROVNANIE
S = [0]*7

#Čiastočné súčty
S[0] = sums[-1]

#Eulerova transformácia (iba postupnosti so striedavými znamienkami)
S[1] = Euler(a_n, N)

#Opakovaný Aitkenov proces
S[2] = sums
while len(S[2])>=3:
    S[2] = aitken(S[2])
S[2] = S[2][-1]

#Opakovaný Aitkenov proces pre podpostupnosť párnych členov
S[3] = sums[0::2]
while len(S[3]) >=3:
    S[3] = aitken(S[3])
S[3] = S[3][-1]

#Opakovaný Aitkenov proces pre podpostupnosť nepárnych členov
S[4] = sums[1::2]
while len(S[4]) >=3:
    S[4] = aitken(S[4])
S[4] = S[4][-1]

#Priemer
S[5] = (S[3] + S[4])/2

#Van Wijngaarden + Euler (iba postupnosti s rovnakými znamienkami)
M_VW = 25
a_VW = VW(a_n, N, M_VW)
S[6] = Euler_(a_VW)

#ODCHÝLKY
deviations = [0]*7
deviations[0] = abs(S[0] - SUM).evalf(20)
deviations[1] = abs(S[1] - SUM).evalf(20)
deviations[2] = abs(S[2] - SUM).evalf(20)
deviations[3] = abs(S[3] - SUM).evalf(20)
deviations[4] = abs(S[4] - SUM).evalf(20)
deviations[5] = abs(S[5] - SUM).evalf(20)
deviations[6] = abs(S[6] - SUM).evalf(20)


#Grafy
#Opakovaný Aitkenov proces
if len(sums)%2 == 0:
    deviations_A = [0]*(int(len(sums)/2))
    Aitkens = [0]*(int(len(sums)/2))
else:
    deviations_A = [0]*(int((len(sums) + 1)/2))
    Aitkens = [0]*(int((len(sums) + 1)/2))

deviations_A[0] = [0]*len(sums)
Aitkens[0] = sums

for i in range(len(sums)):
    deviations_A[0][i] = abs(sums[i] - SUM).evalf(20)

for i in range(1, len(Aitkens)):
    Aitkens[i] = aitken(Aitkens[i - 1])
    deviations_A[i] = [0]*len(Aitkens[i])
    for p in range (len(deviations_A[i])):
        deviations_A[i][p] = abs(Aitkens[i][p] - SUM).evalf(20)
        
plt.figure()
for i in range(len(deviations_A)):
    os_x = np.linspace(len(deviations_A[0]) - len(deviations_A[i]) , len(deviations_A[0]) - 1, len(deviations_A[i]))
    plt.scatter(os_x, deviations_A[i], marker = 'X')
pyplot.yscale('log')
plt.xlabel('n', fontsize = 15)
plt.ylabel('$|S_{n}-L|$', fontsize = 15)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.grid()

#Eulerova transformácia
deviations_Euler = [0]*2
deviations_Euler = [[0]*len(sums) for i in range(2)]
for i in range(len(sums)):
    deviations_Euler[0][i] = abs(sums[i] - SUM).evalf(20)
    deviations_Euler[1][i] = abs(Euler(a_n, i) - SUM).evalf(20)
os_x = np.linspace(0, len(sums) - 1, len(sums))
plt.figure()
plt.scatter(os_x, deviations_Euler[0], marker = 'X')
plt.scatter(os_x, deviations_Euler[1], marker = 'X')
pyplot.yscale('log')
plt.xlabel('n', fontsize = 15)
plt.ylabel('$|S_{n}-L|$', fontsize = 15)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.grid()

#Shanksova transformácia
k = 2
n = 3
def sucty(n):
    return Sum(a_n(x), (x, 0, n)).doit()
postupnosti = shanks(sucty, n, k)
chyby = np.zeros([k + 1, n + 1])
for i in range(k + 1):
    for p in range(n + 1):
        chyby[i][p] = abs(float(postupnosti[i][p]) - SUM).evalf(20)