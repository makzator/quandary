import numpy

delta = 0.001
cf = 2
f = open("fine_eig_decay.txt","r")
real = 0
imag = 0
for line in f:
    words = line.split()
    real = float(words[0])
    imag = float(words[2])
    break
f.close()
lamb1 = real+imag*complex('j')
#lamb1 is now the first eigenvalue of M (complex)

#formula for implicit midpoint rule eigenvalues
lamb = (2+lamb1*delta)/(2-lamb1*delta)


"""
f = open("coarse_eig.txt","r")
real = 0
imag = 0
for line in f:
    words = line.split()
    real = float(words[0])
    imag = float(words[2])
    break
f.close()
#lamb1 = real+imag*complex('j')
"""
#corresponding coarse eigenvalue
mu = (2+lamb1*cf*delta)/(2-lamb1*cf*delta)

#number of coarse grid time steps
nc = 2.0/(delta*cf)


d1 = abs(mu-lamb**cf)

#equation (1) from ben's email
cv = (d1**2/(1.0-abs(mu)**2+3/nc**2))**0.5

print(cv)
