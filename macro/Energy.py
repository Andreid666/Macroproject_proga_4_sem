import numpy as np
import matplotlib.pyplot as plt

f = open('results/Energy3.xyz', 'r')
y=[]
x=[]
try:
    k = 0
    text = f.readlines()
    for i in text:
        #print(i)
        y.append(float(i))
        x.append([k])
        k += 1
finally:
    f.close()

plt.plot(x, y)
plt.show()
