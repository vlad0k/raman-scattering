from math import sqrt, exp, cos
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


i = 1j

# params / input


wave_length = 2.4 * 10**-12

# constants

e = 4.8 * 10**-10 
c = 2.9979 * 10**10
m = 9.1094 * 10**-28
pi = 3.14

beta = 1 # берем 1

gamma0 = 4

uii = c * sqrt(1 - (1 / gamma0)**2)
w = 1 * 10**12 # -s # частота волны
Omega = 10 ** 11 # -s
T = 6 * 10**-12
t0 = 0
t1 = T
n = 2000

kii = (w - Omega) / uii

lmbda = 2 * pi / kii

Eii = 10**-4

# equations

def z_n(z0, t): # z невозмущенное
  return z0 + uii * (t - t0)
  
def E(z,t): # wave field
  # return 0.5 * abs(Eii) * сexp(-i * w * t + i * kii * z) * 2
  return 0.5 * abs(Eii) * cos(-w * t + kii * z) * 2

def W (u): # relativistic kinetic energy u - current speed 
  return (m * c**2) / sqrt(1 - u**2 / c**2) - m * c**2

def system_solve(z0, t_list):
  def F(s,t): # s - vector [z, v]
    dzdt = s[1]
    dvdt = (e / m) * (beta ** 2) * ((1 - (s[1] / c )** 2) ** 1.5) * E(s[0], t) - (Omega ** 2) * (s[0] - z_n(z0,t))
    return [dzdt, dvdt]

  s0=[z0, uii]
  s = odeint(F, s0, t_list)
  return s


def find_electron_energy(z0, T):  #z0 - зависит от положения на волне
  U = [u[1] for u in system_solve(z0, T)]
  return [W(u) for u in U]


# start!!!

Z = np.linspace(0, lmbda, n, endpoint=False) # generate z0 for 1000 electrons
T_list = np.linspace(t0,t1, n, endpoint=False) # 1s interval for 30s 

electron_energies = [find_electron_energy(z0, T_list) for z0 in Z]

W_mean = [np.mean(i) for i in np.array(electron_energies).transpose()]

plt.plot(T_list,W_mean)
plt.show()



# ignore plz
# def acceleration (t): #u - current speed
  # return (e / m) * (beta ** 2) * ((1 - (u(t) ** 2) / (c ** 2)) ** 1.5) * E(t) - (Omega ** 2) * (z(t) - z_n(t))