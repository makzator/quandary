from qutip import*
import numpy as np
import matplotlib.pyplot as plt
import time

def H1_coeff(t, args):
    return args['Jkl'] * (np.cos(args['eta']*t))

def H2_coeff(t, args):
    return args['Jkl'] * (np.sin(args['eta']*t))

nsteps = 1000
tlist = np.linspace(0,0.01,nsteps+1) # time in micro sec (used to be 1001 steps to time 0.01)

Nq = 3
Nc = 3
a = tensor(destroy(Nq), qeye(Nc))
c = tensor(qeye(Nq), destroy(Nc))

H1 = (a.dag()*c + a*c.dag())
H2 = 1j*(a.dag()*c - a*c.dag())
H3 = a + a.dag()
H4 = 1j*(a - a.dag())
H5 = c + c.dag()
H6 = 1j*(c - c.dag())

chi_a = 3.58    # MHz
chi_ac = 10     # MHz
chi_c = 4.26    # MHz
freq_a = 3.0 # fundamental frequency for the qudit, MHz
freq_c = 3.0 # fundamental frequency for the cavity, MHz
rfreq_a = 3.0  # rotation frequency for the qudit, MHz
rfreq_c = 3.0  # rotation frequency for the cavity, MHz
detuning_a = freq_a - rfreq_a # detuning frequency for the qudit, MHz
detuning_c = freq_c - rfreq_c # detuning frequency for the cavity, MHz
r1_a = 1./95.33    # use 95.33us T1 time
r2_a = 1./134.65   # use 134.65us T2 time
r1_c = 1./89.46    # use 89.46us T1 time
r2_c = 1./152.11   # use 152.11us T2 time

# drive from |2> towards 1> state
args = {'Jkl': 2*np.pi*chi_ac, 'eta': 2*np.pi*(rfreq_a - rfreq_c)}
# rho_ini = tensor(basis(Nq,2),basis(Nc,0)) # starting Alice from |2>
axc0 = tensor(basis(Nq,2), basis(Nc,0))
rho_ini = axc0 * axc0.dag()
print("initial state: ", rho_ini)

# Note the factor 2*pi in H
H0 = (detuning_a*a.dag()*a + detuning_c*c.dag()*c - chi_a/2*(a.dag()*a.dag())*(a*a) \
    -chi_c/2*(c.dag()*c.dag())*(c*c))*2*np.pi + H3 + H4 + H5 + H6
H = [H0, [H1, H1_coeff], [H2, H2_coeff]]

c_ops = [np.sqrt(r1_a)*a, np.sqrt(r1_c)*c, np.sqrt(r2_a)*(a.dag()*a), np.sqrt(r2_c)*(c.dag()*c)]
e_ops = [a.dag()*a, c.dag()*c]
#e_ops = [sigmay(), sigmaz()]

options = Options(atol=1e-16, rtol=1e-16)
print(options)

start_time = time.time()
# output = mesolve(H, rho_ini, tlist, [], e_ops, args = args)
output = mesolve(H, rho_ini, tlist, c_ops, [], args = args, options=options)
print("--- %s seconds ---" % (time.time() - start_time))
output.times

print("last real state: ", np.transpose(np.real(output.states[1000])).ravel()) 
print("last imaginary state: ", np.transpose(np.imag(output.states[1000])).ravel()) 


# print(np.real(output.states[1000]*output.states[1000].dag()))
# print(np.imag(output.states[1000]*output.states[1000].dag()))

# # Print the state rho(T) 
# print("State at final time ", tlist[-1])
# rho = output.states[-1]
# # print(np.real(rho))
# print(operator_to_vector(rho))

# # Read in quandary expected energy file
# def getdata(filename):
#   # open file for read only
#   ufile=open(filename,'r')
#   times=[]
#   data=[]
#   # parse the lines of the file
#   for i,line in enumerate(ufile):
#     if (line.startswith('#')):     # ignore lines that begin with '#'
#         continue
#     line = line.split()
#     times.append(float(line[0]))
#     data.append(float(line[1]))
#   #close the file
#   ufile.close()
#   return times, data 

# times, quan_expected0 = getdata("/Users/guenther5/Numerics/quandary/results/ibm_q5q6/data_out/expected0.iinit-001.rank0000.dat")
# times, quan_expected1 = getdata("/Users/guenther5/Numerics/quandary/results/ibm_q5q6/data_out/expected1.iinit-001.rank0000.dat")

# # Plot expected energy level over time
# fig, ax = plt.subplots(figsize=(12,6))
# ax.plot(tlist, output.expect[0], 'r')
# ax.plot(tlist, output.expect[1], 'b')
# ax.plot(times, quan_expected0, 'g')
# ax.plot(times, quan_expected1, 'y')
# ax.legend(("Qutip0", "Qutip1", "Quandary0", "Quandary1"))
# ax.set_xlabel('time')
# ax.set_ylabel('expectation energy level');
# plt.show(block=True)
