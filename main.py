from scipy.constants import e, c, m_e as m
from math import sqrt, exp

i = complex(0, 1)

# params

Eii = 0
Omega = 0
uii = 0
t = 0

# constants

t0 = 0
z0 = 0
beta = 0
kii = 0

# equations

def z(t):
  return z0 + uii * (t - t0)

def w(t): 
  return kii * uii

def E(t): # wave field
  return 0.5 * abs(Eii) * exp( -i * w(t) * t + i * kii * z(t))

def W (): # relativistic kinetic energy
  return (m * c**2) / sqrt(1 - uii**2 / c**2)

def acceleration ():
  return (e / m) * (beta ** 2) * ((1 - (uii ** 2) / (c ** 2)) ** 1.5) * E(Eii) - (Omega ** 2) * (z(t) - z0)
