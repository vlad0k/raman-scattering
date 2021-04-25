from math import sqrt, exp
import numpy as np
from scipy.integrate import odeint

i = complex(0, 1)



# params / input

Eii = 0
Omega = 0
uii = 1 #начальная скорость
t = 0

wave_length = 2.4 * 10**-12

# constants

e = 4.8 * 10**-10 
c = 2.9979 * 10**10
m = 9.1094 * 10**-28

t0 = 0
beta = 1 
kii = 0 
w = 0 # частота волны



# equations

def z_n(z0, t): # z невозмущенное
  return z0 + uii * (t - t0)
  

def E(z,t): # wave field
  return 0.5 * abs(Eii) * exp( -i * w * t + i * kii * z)


def W (u): # relativistic kinetic energy u - current speed 
  return (m * c**2) / sqrt(1 - u**2 / c**2) - m * c**2


def system_solve(z0, t_list):
  def F(s,t): # s - vector [z, v]
    print(t, s[0], s[1])
    dzdt = s[0]
    dvdt = (e / m) * (beta ** 2) * ((1 - (s[1] ** 2) / (c ** 2)) ** 1.5) - (Omega ** 2) * (s[0] - z_n(z0,t))

    # cant hanlde complex value of E(z,t)
    # dvdt = (e / m) * (beta ** 2) * ((1 - (s[1] ** 2) / (c ** 2)) ** 1.5) * E(s[0], t) - (Omega ** 2) * (s[0] - z_n(t))
    
    return [dzdt, dvdt]

  s0=[z0, uii]
  print(t_list)
  s = odeint(F, s0, t_list)
  return s


def find_electron_energy(z0 = 0, T = np.linspace(0,30, 30)):  #z0 - зависит от положения на волне
  U = [u[1] for u in system_solve(z0, T)]
  W_list = []
  for u in U:
    W_list.append(W(u))

  return W_list


# start!!!

Z0 = np.linspace(0, wave_length, 1000, endpoint=False) # generate z0 for 1000 electrons
T = np.linspace(0,30, 30, endpoint=False) # 1s interval for 30s 

for z0 in Z0:
  print(find_electron_energy(z0, T)) # prints all electron energies for every electron





# ignore plz
# def acceleration (t): #u - current speed
  # return (e / m) * (beta ** 2) * ((1 - (u(t) ** 2) / (c ** 2)) ** 1.5) * E(t) - (Omega ** 2) * (z(t) - z_n(t))
