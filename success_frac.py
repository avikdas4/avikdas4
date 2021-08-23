from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
import glob
from scipy.interpolate import interp1d
import numpy as np
import sys

###############

file_no = input("Enter the number of simulated light curves in your directory:")

###############

file = open('gamma_whole_sort.txt') # Your file name!
lines = file.readlines()
data = np.loadtxt(lines)

Time = data[:,0]
flux = data[:,1]


start_time = np.min(Time)
end_time  = np.max(Time)

#print(start_time, end_time)
#sys.exit()

flux1 = interp1d(Time, flux, fill_value="extrapolate")

reg_Time = np.arange(start_time, end_time, 10.0)
reg_flux = flux1(reg_Time)

plt.plot(Time, flux,'bo')
plt.plot(reg_Time, reg_flux, 'm^')

#print(reg_flux.size)
#print(flux.size)

np.savetxt("gamma_resample.dat", np.transpose([reg_Time, reg_flux]), fmt="%e\t%e") # you can change the file name here!

plt.show()
###################

Time_r = reg_Time 
flux_r = reg_flux

###################

f_obs, pxx_den_obs = signal.periodogram(flux_r,fs=1,scaling='spectrum') #,detrend='linear')
plt.loglog(f_obs, pxx_den_obs,'bo')

freq_number = f_obs.size - 1

f_obs_new = np.zeros(freq_number)
psd_obs = np.zeros(freq_number)
for j in range(freq_number):
#print(f[1])
	f_obs_new[j] = f_obs[j+1]
	psd_obs[j] = pxx_den_obs[j+1]


###################

t = np.zeros((file_no,freq_number,2)) #20
t_sq = np.zeros((file_no,freq_number,2))

for i in range(file_no):

	file1 = open("lightcurve_beta10_%d.txt" % i)
	lines = file1.readlines()
	data1 = np.loadtxt(lines)

	T = data1[:,0]
	flux = data1[:,1]
#	fluxerr= data[:,2]*(10**(6))



	f, pxx_den = signal.periodogram(flux,fs=1,scaling='spectrum')#,detrend='linear')
	plt.loglog(f, pxx_den,'bo')
	f, pxx_den = signal.periodogram(flux,fs=1,scaling='spectrum')#,detrend='linear')
#plt.loglog(f, pxx_den,'bo')

	f_new = np.zeros(freq_number)
	pxx_den_new = np.zeros(freq_number)
	data = np.zeros((freq_number, freq_number))

	for j in range(freq_number):
#print(f[1])
		f_new[j] = f[j+1]
		pxx_den_new[j] = pxx_den[j+1]


	t[i,:,0] = f_new[:]
	t[i,:,1] = pxx_den_new[:]
#	t_sq[i,:,1] = data1[:,1]*data1[:,1]

#print(t.shape)
#sys.exit()
mean_P = np.mean(t,axis=0)
std_P = np.std(t,axis=0)


chi_sim = np.zeros(file_no)
chi_obs = np.zeros(file_no)
for i in range(file_no):
	C = (((t[i,:,1] - mean_P[:,1])**2)/(std_P[:,1])**2)
	C_sum = np.sum(C)	
	chi_sim[i] = C_sum
	K = (((psd_obs - mean_P[:,1])**2)/(std_P[:,1])**2)
	K_sum = np.sum(K)
	chi_obs[i] = K_sum


A = chi_sim/chi_obs

r = 1
count = 0

for j in A:
	if j >= r:
		count = count + 1

success_frac = (float(count) / chi_sim.size)*100

print("success fraction in percentage:" + str(success_frac))








