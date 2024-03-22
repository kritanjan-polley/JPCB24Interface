import random
import math

mu = 0.0
rt = 8.314*300
m = 3*16*1e-3
std = math.sqrt(rt/m)*1e-5

a = -abs(random.gauss(mu,std))
print(a)
