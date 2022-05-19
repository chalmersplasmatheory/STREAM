import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d


hf = h5py.File('65108.h5', 'r')
print(hf.keys())
hf = h5py.File('firdens.h5', 'r')
print(hf.keys())
'''
t_FIR = np.array(hf.get('FIR').get('x'))
FIR = np.array(hf.get('FIR').get('z'))
#i0 = np.argmin(np.abs(t_FIR+0.02))
plt.plot(t_FIR, FIR)# - FIR[i0])
plt.ylabel('Line-integrated (FIR)')
plt.xlabel('t')
plt.xlim([-0.02, 0.1])
plt.show()

'''
hf.close()
hf = h5py.File('TCV65108_old.h5', 'r')
print(hf.keys())
'''
t_CIII = np.array(hf.get('C-III').get('x'))
cIII = np.array(hf.get('C-III').get('z'))
plt.plot(t_CIII, cIII)
plt.ylabel('C-III')
plt.xlabel('t')
plt.xlim([-0.02, 0.05])
plt.show()

t_H = np.array(hf.get('H-alpha').get('x'))
Halpha = np.array(hf.get('H-alpha').get('z'))
plt.plot(t_H, Halpha)
plt.ylabel('H-alpha')
plt.xlabel('t')
plt.xlim([-0.02, 0.05])
plt.show()

t_FIR = np.array(hf.get('Line-integrated (FIR)').get('x'))
FIR = np.array(hf.get('Line-integrated (FIR)').get('z'))
i0 = np.argmin(np.abs(t_FIR+0.02))
plt.plot(t_FIR, FIR - FIR[i0])
plt.ylabel('Line-integrated (FIR)')
plt.xlabel('t')
plt.xlim([-0.02, 0.1])
plt.show()


t_loop = np.array(hf.get('Loop voltage').get('x'))
V_loop = np.array(hf.get('Loop voltage').get('z'))
V_loop_s = savgol_filter(V_loop, 105, 3)
plt.plot(t_loop, -V_loop, t_loop, -V_loop_s)
plt.xlim([-0.01, 0.39])
plt.ylabel('Loop voltage')
plt.xlabel('t')
plt.ylim([-1, 10])
plt.show()

t_maxThomson = np.array(hf.get('Maximum (Thomson)').get('x'))
maxThomson = np.array(hf.get('Maximum (Thomson)').get('z'))
plt.plot(t_maxThomson, maxThomson)
plt.ylabel('Maximum (Thomson)')
plt.xlabel('t')
plt.xlim([-0.02, 0.05])
plt.show()
'''
t_fluxD = np.array(hf.get('Particle flux (D2)').get('x'))
fluxD = np.array(hf.get('Particle flux (D2)').get('z'))
plt.plot(t_fluxD, fluxD)
plt.ylabel('Particle flux (D2)')
plt.xlabel('t')
plt.xlim([-0.01, 0.4])
plt.show()

'''
t_Ip = np.array(hf.get('Plasma current').get('x'))
Ip = np.array(hf.get('Plasma current').get('z'))
plt.plot(t_Ip, Ip)
plt.ylabel('Plasma current')
plt.xlabel('t')
plt.xlim([-0.02, 0.05])
plt.show()

t_Vp = np.array(hf.get('Plasma volume').get('x'))
V_p_osc = np.array(hf.get('Plasma volume').get('z'))
V_p = savgol_filter(V_p_osc, 45, 3)
a_vec = np.sqrt(V_p / (2 * np.pi ** 2 * 0.89))
V_vessel = 4.632
V_initfun = interp1d(np.append(np.array([-0.02]), t_Vp), np.append(np.array([V_vessel]), V_p), 'cubic')
t_Vp2 = np.linspace(-0.02, t_Vp[-1])
V_p2 = V_initfun(t_Vp2)
plt.plot(t_Vp2, V_p2, t_Vp, V_p, t_Vp, V_p_osc)
plt.ylabel('Minor radius')
plt.xlabel('t')
#plt.xlim([-0.02, 0.05])
plt.show()

t_B = np.array(hf.get('Toroidal magnetic field').get('x'))
B = np.array(hf.get('Toroidal magnetic field').get('z'))
plt.plot(t_B, B)
plt.ylabel('Toroidal magnetic field')
plt.xlabel('t')
#plt.xlim([-0.02, 0.05])
plt.show()
'''
hf.close()