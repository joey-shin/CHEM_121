# import packages for basic math, plotting, linear algebra, etc.
from numpy import *
from numpy.linalg import *
from numpy.random import *
from matplotlib.pyplot import *
from scipy.special import binom, erf, erfc

def s(uA,uB,alpha):
    return sqrt(0.5 * pi / alpha) * exp(-0.5 * alpha * (uA - uB) ** 2)

def f(uA,uB,alpha):
    return 0.5*sqrt(0.5 * pi / alpha) * exp(-0.5 * alpha * (uA - uB) ** 2) \
           * (alpha - alpha**2 * (uA - uB)**2 + \
                             0.25*(1/alpha + (uA + uB)**2 ))

def g(uA,uB,alpha):
    return sqrt(0.5 * pi / alpha) * exp(-0.5 * alpha * (uA - uB) ** 2) \
           * (3 / (16 * alpha ** 2) + \
               (3 / (8 * alpha)) * (uA + uB) ** 2 \
               + (1 / 16) * (uA + uB) ** 4)
