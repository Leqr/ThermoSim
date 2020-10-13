import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

#loading the solutions
ffluid = open("../build/debug/fluid.txt")

solsf = []

for line in ffluid:
    a = line.strip().split(',')
    af = []
    for val in a:
        af.append(float(val))
    solsf.append(af)

fsolid = open("../build/debug/solid.txt")

solss = []

for line in fsolid:
    a = line.strip().split(',')
    af = []
    for val in a:
        af.append(float(val))
    solss.append(af)


#animation of the solution
nCells = 100
nTimeStep = 1000

fig = plt.figure()
ax = plt.axes(xlim=(0, nCells),ylim = (0,900))
line, = ax.plot([], [], lw=2)

def ini():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(1, nCells, nCells)
    y = solss[i]
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=ini,
                               frames=nTimeStep, interval=5,blit = True)
plt.show()
