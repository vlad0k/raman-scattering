from math import sqrt, exp, cos, pi
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from cmath import exp as cexp


i = 1j

# params / input


wave_length = 2.4 * 10**-12

# constants

e = -4.8 * 10**-10 
c = 2.9979 * 10**10
m = 9.1094 * 10**-28

beta = 1 # берем 1

gamma0 = 4

uii = c * sqrt(1 - (1 / gamma0)**2)
w = 1 * 10**12 # -s # частота волны
Omega = 10 ** 11 # -s
delta_w = Omega / 5
T = 2 * pi / w

t0 = 0
t1 = 10*T
n = 1000

kii = (w + Omega) / uii

lmbda = 2 * pi / kii

Eii = 10**-4

# equations

def fi(z): 
  kii * z

def z_n(z0, t): # z невозмущенное
  return z0 + uii * (t - t0)
  
def E(z,t): # wave field
  return 0.5 * abs(Eii) * cos(-w * t + kii * z) * 2

def W (u): # relativistic kinetic energy u - current speed 
  return (m * c**2) / sqrt(1 - u**2 / c**2) - m * c**2


def z_wawe(z):
  return (( ((e / m) * beta * gamma0**-3 * 0.5 * Eii) /( (-i * w + delta_w + kii * uii)**2 + Omega**2)) * 2* cexp(i*kii*z)).real * 2


def u_wawe(z):
  return (( ((e / m) * beta * gamma0**-3 * 0.5 * Eii * (-i * w + delta_w - i * kii * uii)) /( (-i * w + delta_w + kii * uii)**2 + Omega**2)) * 2* cexp(i*kii*z)).real * 2


def system_solve(z0, t_list):
  def F(s,t): # s - vector [z, v]
    dzdt = s[1]
    dvdt = (e / m) * (beta ** 2) * ((1 - (s[1] / c )** 2) ** 1.5) * E(s[0], t) - (Omega ** 2) * (s[0] - z_n(z0,t))
    return [dzdt, dvdt]

  s0=[z_wawe(z0), u_wawe(z0)]
  s = odeint(F, s0, t_list, tcrit=4)
  return s

def find_electron_energy(z0, T):  #z0 - зависит от положения на волне
  U = [u[1] for u in system_solve(z0, T)]
  return [W(u) for u in U]



# start!!!


Z = np.linspace(0, lmbda, n, endpoint=False) # generate z0 for 1000 electrons
T_list = np.linspace(t0,t1, n, endpoint=False) # 1s interval for 30s 

electron_energies = [find_electron_energy(z0, T_list) for z0 in Z]

W_mean = [np.mean(i) for i in np.array(electron_energies).transpose()] #энергия пучка

E0= [] # начальная энергия для каждого электрона для каждого момента времени
for t in T_list:
  index = len(E0)
  E0.append([])
  for z in Z:
    E0[index].append(E(z, t))

E = [np.mean(i) for i in np.array(E0).transpose()] # начальная энергия поля в каждый момент времени


y = [-(E[i] - W_mean[i]) for i in range(len(W_mean))] # изменение энергии поля, где (0.5 * e * E[i]) - начальная энергия, 

plt.plot(T_list,y)
plt.show()



# ignore plz
# def acceleration (t): #u - current speed
  # return (e / m) * (beta ** 2) * ((1 - (u(t) ** 2) / (c ** 2)) ** 1.5) * E(t) - (Omega ** 2) * (z(t) - z_n(t))