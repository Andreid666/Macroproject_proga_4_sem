import numpy as np
import matplotlib.pyplot as plt


def plot_c(c, x_axis, y_axis, color, v_x):

    #x_axis = []
    t = max(v_x) - min(v_x)
    for i in range(0, c+1):
        x_axis.append(min(v_x) + t/(2*c) + i * t/c)
        y_axis.append(0)

    for k in v_x:
        l = int(c*((k - min(v_x))/t))
        #print(l)
        for j in range (l - 5, l + 5):
            if (j >= 0 and j < c + 1):
                y_axis[j] += 1

    new_x = []
    new_y =[]

    for i in range(0, len(y_axis)):
        if y_axis[i] != 0 and x_axis[i] != 0:
            new_x.append(x_axis[i])
            new_y.append(y_axis[i])
    #print(new_y)
    #print(new_x)
    y_axis = new_y
    #print(y_axis)
    x_axis = new_x
    #print(x_axis)

    return new_x, new_y

    #plt.plot(x_axis, y_axis, color)

def calcs(o, f):
    v_x = []
    text = f.readlines()
    #o = 390
    size = int(text[0])
    for i in range (1 + o * (size+1), 1 + o* (size+1) + size):
        v_x.append(float(text[i]))
    return v_x



#v_x = []

x_axis = []
y_axis = []
#x_axis1 = []
#x_axis2 = []

c_0 = 1000
c_1 = 40
c_2 = 80

f = open('results/Maxwell2.xyz', 'r')
try:
    v_x = calcs(50, f)

finally:
    #for i in v_x:
    #print(i)
    #print("asda")
    f.close()




number_of_particles = 125000
p = plot_c(c_0, x_axis, y_axis, 'r', v_x)
#print(x_axis)
#print(y_axis)
'''
x = np.empty(len(p[0]))
y = np.empty(len(p[1]))
for i in range (0, len(p[0])):
    if (p[0][i] == 0 or p[1][i] == 0):
        print ("WHYTHEFUCK")
    x[i] = p[0][i]

    y[i] = p[1][i]
'''
x = np.array(p[0])
y = np.array(p[1])


#t = np.linspace(-1.2*1e-14, 0, len(x))

#w = np.power(x, 2)
z = 0.5*np.sum(x)/number_of_particles
#w = w
print(z)

plt.plot(x,  (np.log(y[0])- np.log(x[0]) + 0*np.log(4*number_of_particles/(np.sqrt(np.pi*np.power(number_of_particles*z, 3))))) - (x-x[0])/z, 'b')

plt.plot(x, np.log(y)- np.log(x), 'r')
#plot_c(c_1, x_axis, y_axis, 'g', v_x)
#plot_c(c_2, x_axis, y_axis, 'b', v_x)

plt.show()