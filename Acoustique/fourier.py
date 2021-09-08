 
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:54:33 2018

@author: maugeais
"""

import numpy as np
import matplotlib.pyplot as plt
# Calcul les coefficients de fourier de 0 à N
# f est une fonction L périodique
# Intégration numérique par la méthode des rectangles
def coeff(f, N, T = 2*np.pi) :
    a = []
    b = []
    
    dx = T/1000
    x = np.arange(0, T, dx)
    print(T)
    
    for n in range(N) :
        a.append(np.sum(f(x+dx/2)*np.cos(2*np.pi*n*(x+dx/2)/T))*dx*2/T)
        b.append(np.sum(f(x+dx/2)*np.sin(2*np.pi*n*(x+dx/2)/T))*dx*2/T)

    return(a, b)

# Construit la série de fourier à partir des coefficients évalué sur x
def series(a, b, x, T) :
    S = np.zeros(len(x))+a[0]/2
  
    for n in range(1, len(a)) :
        S += a[n]*np.cos(2*np.pi*n*x/T)+b[n]*np.sin(2*np.pi*n*x/T)
        
    return(S)


if __name__ == '__main__':
    
    f1 = lambda x: np.mod(x, 2*np.pi) > np.pi #np.cos(3*x)
    f2 = lambda x: np.mod(x, 2*np.pi)
    f4 = lambda x: np.mod(x, 2*np.pi)*(2*np.pi-np.mod(x, 2*np.pi))
    
    g2 = lambda x: abs(np.mod(x, 2*np.pi)-np.pi)
    
    f = f4
    
    plt.ion()
    # Calcul et affichage des coefficients de fourier
    a, b = coeff(f, 100)
    plt.figure(0)
    plt.plot(a, '-ro', linewidth=0)
    plt.plot(b, '-bo', linewidth=0)
    plt.show()
    
    # Calcul et affichage de la fonction et de sa série de Fourier
    dx = 2*np.pi/100
    x = np.arange(-2*np.pi, 4*np.pi, dx)
    
    S = series(a, b, x)
    plt.figure(1)
    plt.plot(x, f(x))
    plt.plot(x, S)
    plt.show()
    
    #input()
            
            
            