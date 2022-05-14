# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 19:28:58 2022

@author: vincz
"""

from sympy import *

init_session(use_latex=True)


#Aitkenov proces
def aitken(A):
    B = [0] * (len(A) - 2)
    for i in range(1, len(A) - 1):
        #B[i-1] = simplify((A[i-1] * A[i+1] - A[i]**2) / (A[i-1] + A[i+1] - 2*A[i]))
        B[i - 1] = (A[i - 1] * A[i + 1] - A[i]**2) / (A[i - 1] + A[i + 1] - 2*A[i])
    return B


#Shanksova transformácia
def shanks(a_n, n, k):
    matice = [0]*(k + 1)
    pole = [0]*(k + 1)
    for i in range(k + 1):
        matice[i] = [0]*(n + 1)
        pole[i] = [0]*(n + 1)
        for p in range(n + 1):
            matice[i][p] = [0]*2
            matice[i][p][0] = [0]*(i + 1)
            matice[i][p][1] = [0]*(i + 1)
            for t in range(i + 1):
                matice[i][p][0][t] = [0]*(i + 1)
                matice[i][p][1][t] = [0]*(i + 1)
            for t in range(i + 1):
                matice[i][p][0][0][t] = a_n(p + t)
                matice[i][p][1][0][t] = 1
            for t in range(1, i + 1):
                for q in range(i + 1):
                    matice[i][p][0][t][q] = a_n(p + t + q) - a_n(p + t + q - 1)
                    matice[i][p][1][t][q] = a_n(p + t + q) - a_n(p + t + q - 1)
            matice[i][p][0] = np.array(matice[i][p][0], dtype = float)
            matice[i][p][1] = np.array(matice[i][p][1], dtype = float)
            pole[i][p] = (np.linalg.det(matice[i][p][0]))/(np.linalg.det(matice[i][p][1]))
    return pole


#Eulerova transformácia 
#(daný predpis sčítavanej postupnosti, výstupom je N-tý súčet transformovanej)
def Euler(a_n, N):
    a = [0]*N
    for i in range(N):
        a[i] = a_n(i)
        
    for_dif_op = [0]*N
    for i in range(N):
        for p in range(i + 1):
            for_dif_op[i] = for_dif_op[i]+((-1)**p)*binomial(i,p)*  abs(a[i - p])
    
    Euler_sum = 0
    for i in range(N):
        Euler_sum = Euler_sum + ((-1)**i)*for_dif_op[i]/(2**(i + 1))
    return Euler_sum

#(daných niekoľko prvých členov postupnosti, výstupom je N-tý súčet transformovanej)
def Euler_(a):
    for_dif_op = [0]*len(a)
    for i in range(len(a)):
        for p in range(i+1):
            for_dif_op[i] = for_dif_op[i]+((-1)**p)*binomial(i,p)*abs(a[i-p])
        
    Euler_sum = 0
    for i in range(len(a)):
        Euler_sum = Euler_sum+((-1)**i)*for_dif_op[i]/(2**(i+1))
    return Euler_sum


#Van Wijngaardenova transformácia
def VW (a_n, N, M):
    A = [0]*N
    b = [0]*N
    for i in range(N):
        b[i] = [0]*M
        for p in range(M):
            b[i][p] = (2**p)*a_n((2**p)*(i + 1) - 1)
        A[i] = ((-1)**i)*sum(b[i])
    return A