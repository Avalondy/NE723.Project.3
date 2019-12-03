#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
from matplotlib.pyplot import figure

filename = "data.for.plotCharacteristicMesh.txt"

#x = np.linspace(0.,20.16,num=225)
#filename = "results/results.test.B.txt"

graph_data = open(filename, 'r').read()
lines = graph_data.split('\n')

x_plot = np.array([])
y_plot = np.array([])

for line in lines:
    if len(line)>0:
        x, y = line.split()
        x_plot = np.append(x_plot, float(x))
        y_plot = np.append(y_plot, float(y))
    elif x_plot.size !=0:
        plt.plot(x_plot, y_plot, marker='.', color="black", markerfacecolor="r")
        x_plot = np.array([])
        y_plot = np.array([])



# plt.scatter(x_plot, y_plot)
# plt.show()
plt.savefig('plotCharacteristicMesh.pdf')
