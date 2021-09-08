#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 08:31:40 2021

@author: maugeais
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile

from scipy.signal import fftconvolve

import fourier

""" Paramètres de l'instrument """
F1 = 72
L = 0.65
c =  1/(2*L*F1)  # = sqrt(T/M))
sigma = L**2*c**2/np.pi**2
eta= 3e-6

""" Conditions initiales : décompostion de Fourier """
N = np.arange(1, 100, dtype=int)

p = 0.1 # triangle de point le plus haut en p
def f(x) :
    y = np.mod(x, 2*L)
    y = y*(y < L)-(2*L-y)*(y >= L)
    z = (y >= 0) * (y/p*(y < p) + (L-y)/(L-p)*(y >= p))
    z -= z[::-1]

    return(z)

a, b = fourier.coeff(f, N[-1]+1, 2*L)

x = np.arange(-L, L, 1e-3)
xp = np.arange(-L, L, 1e-3)
S = fourier.series(a, b, xp, 2*L)
plt.figure(1)
plt.plot(xp, f(xp))
plt.plot(xp, S, label='N = '+str(N[-1]))
plt.show()
plt.legend()

plt.figure(2)
plt.plot(N*b[1:])

""" Construction du signal au chevalet, donné par la dérivée en première en espace
parce que c'est la pente de la déformation qui compte !!! """
omega = (1j*eta*N**2+np.sqrt(4*sigma*N**2-eta**2*N**4))/(-2*sigma)

T = 10
t = np.arange(0, T, 1/44100)
w = np.zeros_like(t)

for n in range(len(N)):
    w += (n+1)*b[n+1]*np.cos(np.real(omega[n])*t)*np.exp(np.imag(omega[n])*t)
    
plt.figure(3)    
plt.plot(t, w)

wavfile.write('signal.wav', 44100, w/max(abs(w)))
plt.figure('spectre')
s = np.fft.fft(w)
f = np.fft.fftfreq(len(w), 1/44100)
I = np.where(np.logical_and(f > 0, f < 2000))[0]
plt.semilogy(f[I], abs(s[I]), label='N = '+str(N[-1]))
plt.legend()

    
""" convolution avec la FRF """
plt.figure('spectreFRF')
rate, FRF = wavfile.read('FRF1.wav')
sFRF = np.fft.fft(FRF[:, 0])
f = np.fft.fftfreq(len(sFRF), 1/rate)
I = np.where(np.logical_and(f > 0, f < 2000))[0]
plt.semilogy(f[I], abs(sFRF[I]))




plt.figure('signalComplet')
w = fftconvolve(w, FRF[:, 0], mode='full')

t = np.arange(0, len(w)/44100, 1/44100)
plt.plot(t, w)

wavfile.write('signalFRF.wav', 44100, w/max(abs(w)))

plt.figure('spectre Complet')
s = np.fft.fft(w)
f = np.fft.fftfreq(len(w), 1/rate)
I = np.where(np.logical_and(f > 0, f < 2000))[0]
plt.semilogy(f[I], abs(s[I]), label='N = '+str(N[-1]))
plt.legend()
