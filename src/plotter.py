import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

def plotMov(nCells,nTimeStep,sol):
    #animation of the solution

    fig = plt.figure()
    ax = plt.axes(xlim=(0, nCells),ylim = (-5,5))
    line, = ax.plot([], [], lw=2)

    def ini():
        line.set_data([], [])
        return line,

    # animation function.  This is called sequentially
    def animate(i):
        x = np.linspace(1, nCells, nCells)
        y = sol[i]
        line.set_data(x, y)
        return line,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=ini,
                                   frames=nTimeStep, interval=5,blit = True)
    plt.show()

def MMS():
    error = np.absolute([np.cos(k*val) for val in np.linspace(0,10,100)]-solss[-1])
    plt.plot(solss[-1])
    plt.show()

if __name__ == "__main__":

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

sol = solss
#plotMov(32,len(sol),sol)
plt.plot(sol[-1])
#plt.plot([np.cos(2*3.14*4/10*x) for x in np.linspace(0,10,8)],label = "cos")
plt.legend()
plt.show()
