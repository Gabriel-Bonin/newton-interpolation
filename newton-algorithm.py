# -*- coding: utf-8 -*-

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

# list of values to be interpolated (x,y)
# below is just an example
xi = np.array([-54, -34, -14, 2, 22, 34])
fi = np.array([102.2, 105.7, 109.2, 109.9, 104, 100])

# divided difference table
title = ['i','xi','fi']
n = len(xi)
ki = np.arange(0,n,1)
table = np.concatenate(([ki],[xi],[fi]),axis=0)
table = np.transpose(table)

# doing divided difference
dfinite = np.zeros(shape=(n,n),dtype=float)
table = np.concatenate((table,dfinite), axis=1)

# generating table
[n,m] = np.shape(table)
diagonal = n-1
j = 3
while (j < m):
    # title for each column
    title.append('F['+str(j-2)+']')

    # every row of each column
    i = 0
    pace = j-2 # starts with 1
    while (i < diagonal):
        denominator = (xi[i+pace]-xi[i])
        numerator = table[i+1,j-1]-table[i,j-1]
        table[i,j] = numerator/denominator
        i = i+1
    diagonal = diagonal - 1
    j = j+1

# divided difference polynomium
dDivided = table[0,3:]
n = len(dfinite)

# rewriting the polynomium with scipy
x = sym.Symbol('x')
polynomium = fi[0]
for j in range(1,n,1):
    factor = dDivided[j-1]
    end = 1
    for k in range(0,j,1):
        end = end*(x-xi[k])
    polynomium = polynomium + end*factor

# simplification by multiplying around (x-xi)
polysimple = polynomium.expand()

# polynomium given numericaly
px = sym.lambdify(x,polysimple)

# defining scale
sample = 101
a = np.min(xi)
b = np.max(xi)
pxi = np.linspace(a,b,sample)
pfi = px(pxi)

# outputs
np.set_printoptions(precision = 4)
print('Divided Difference Table')
print([title])
print(table)
print('dDivided: ')
print(dDivided)
print('polynomium: ')
print(polynomium)
print('simplified polynomium: ' )
print(polysimple)

# Graphic
plt.plot(xi,fi,'o', label = 'Points')
plt.plot(pxi,pfi, label = 'Polynomium')
plt.legend()
plt.xlabel('xi')
plt.ylabel('fi')
plt.title('Differences Division - Newton')
plt.show()
