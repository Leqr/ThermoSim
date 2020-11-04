import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from matplotlib.animation import FuncAnimation

def plotMov(nCells,nTimeStep,sol):
    #animation of the solution

    fig = plt.figure()
    ax = plt.axes(xlim=(0, nCells-1),ylim = (-2,2))
    line, = ax.plot([], [], lw=2)

    def ini():
        line.set_data([], [])
        return line,

    # animation function.  This is called sequentially
    def animate(i):
        x = np.linspace(0, nCells-1, nCells)
        y = sol[i]
        line.set_data(x, y)
        return line,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=ini,
                                   frames=nTimeStep, interval=5,blit = True)
    plt.show()

def plotMovDouble(nCells,nTimeStep,sol):
    #animation of the solution
    sol_fluid = sol[0]
    sol_solid = sol[1]

    fig = plt.figure()
    ax = plt.axes(xlim=(0, nCells-1),ylim = (100,1000))
    x, ydata0, ydata1 = [], [], []
    ln0, = plt.plot([], [], 'r',label = "Fluid Temperature")
    ln1, = plt.plot([], [], 'b',label = "Solid Temperature")

    def init():
        ln0.set_data(x,ydata0)
        ln1.set_data(x,ydata1)
        return ln0, ln1

    def animate(i):
        x = np.linspace(0, nCells-1, nCells)
        ydata0 = sol_fluid[i]
        ydata1 = sol_solid[i]
        ln0.set_data(x, ydata0)
        ln1.set_data(x, ydata1)
        return ln0, ln1,

    ani = FuncAnimation(fig, animate, init_func=init,
                                   frames=nTimeStep, interval=5, blit = True)
    plt.legend()
    plt.show()

def OVSNC():
    #ovs = open("../Data/OVSNonCoupled3/ovsNCn=1,Pe=1.txt")
    ovs = open("../build/debug/ovs.txt")

    ovsdata = []

    for line in ovs:
        a = line.strip().split(',')
        af = []
        for val in a:
            af.append(float(val))
        ovsdata.append(af)

    fig, axs = plt.subplots(2, 2)
    vals = [8,16,32,64,128,256,512,1024]
    '''
    axs[0, 0].plot([np.log(10/n) for n in vals],[np.log(val) for val in ovsdata[2]], 'tab:green',marker=".")
    axs[0, 0].set_title('L1 norm for Ts, n = 4')
    axs[0, 0].set_ylabel('Log(E)')
    axs[0, 0].set_xlabel('Log(h)')
    axs[0, 0].set_ylim([-8,10])


    axs[0, 1].plot([np.log(10/n) for n in vals],[np.log(val) for val in ovsdata[3]], 'tab:orange',marker=".")
    axs[0, 1].set_title('Linf norm for Ts, n = 4')
    axs[0, 1].set_ylabel('Log(E)')
    axs[0, 1].set_xlabel('Log(h)')
    axs[0, 1].set_ylim([-8,10])

    #y = [np.log(np.abs((ovsdata[2][i+1]-ovsdata[2][i])/(ovsdata[2][i]-ovsdata[2][i-1])))/np.log(1/2) for i in range(1,len(ovsdata[2])-1)]
    y = [np.log(ovsdata[2][i]/ovsdata[2][i-1])/np.log(1/2) for i in range(1,len(ovsdata[2]))]

    axs[1, 0].plot([np.log(10/n) for n in vals[1:]],y ,'tab:green',marker=".")
    axs[1, 0].set_ylabel('p')
    axs[1, 0].set_xlabel('Log(h)')
    axs[1, 0].set_ylim([0,4])

    #y2 = [np.log(np.abs((ovsdata[3][i+1]-ovsdata[3][i])/(ovsdata[3][i]-ovsdata[3][i-1])))/np.log(1/2) for i in range(1,len(ovsdata[3])-1)]
    y2 = [np.log(ovsdata[3][i]/ovsdata[3][i-1])/np.log(1/2) for i in range(1,len(ovsdata[3]))]


    axs[1, 1].plot([np.log(10/n) for n in vals[1:]],y2, 'tab:orange',marker=".")
    axs[1, 1].set_ylabel('p')
    axs[1, 1].set_xlabel('Log(h)')
    axs[1, 1].set_ylim([0,4])


    '''
    fig, axs = plt.subplots(2, 2)

    axs[0, 0].plot([np.log(10/n) for n in vals],[np.log(val) for val in ovsdata[0]], 'tab:green',marker=".")
    axs[0, 0].set_title('L1 norm for Tf, n = 1, Pe = 1')
    axs[0, 0].set_ylabel('Log(E)')
    axs[0, 0].set_xlabel('Log(h)')
    axs[0, 0].set_ylim([-8,10])

    axs[0, 1].plot([np.log(10/n) for n in vals],[np.log(val) for val in ovsdata[1]], 'tab:orange',marker=".")
    axs[0, 1].set_title('Linf norm for Tf, n = 1, Pe = 1')
    axs[0, 1].set_ylabel('Log(E)')
    axs[0, 1].set_xlabel('Log(h)')
    axs[0, 1].set_ylim([-8,10])

    #y = [np.log(np.abs((ovsdata[0][i+1]-ovsdata[0][i])/(ovsdata[0][i]-ovsdata[0][i-1])))/np.log(1/2) for i in range(1,len(ovsdata[0])-1)]
    y = [np.log(ovsdata[0][i]/ovsdata[0][i-1])/np.log(1/2) for i in range(1,len(ovsdata[0]))]
    axs[1, 0].plot([np.log(10/n) for n in vals[1:]],y ,'tab:green',marker=".")
    axs[1, 0].set_ylabel('p')
    axs[1, 0].set_xlabel('Log(h)')
    axs[1, 0].set_ylim([0,4])

    #y2 = [np.log(np.abs((ovsdata[1][i+1]-ovsdata[1][i])/(ovsdata[1][i]-ovsdata[1][i-1])))/np.log(1/2) for i in range(1,len(ovsdata[1])-1)]
    y2 = [np.log(ovsdata[1][i]/ovsdata[1][i-1])/np.log(1/2) for i in range(1,len(ovsdata[1]))]

    axs[1, 1].plot([np.log(10/n) for n in vals[1:]],y2, 'tab:orange',marker=".")
    axs[1, 1].set_ylabel('p')
    axs[1, 1].set_xlabel('Log(h)')
    axs[1, 1].set_ylim([0,4])

    plt.legend()
    plt.show()

def plotvscos(sol,ncells,n):
    plt.plot(sol[-1])
    plt.plot([np.cos(2*3.14*n/10*x) for x in np.linspace(0,10,ncells)],label = "real cos")
    plt.suptitle("nC=32,Ts,n=1,ds=5.1e-5")
    plt.legend()
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

sol = solsf

plotMovDouble(256,len(solss),[solsf,solss])
#plotMov(512,len(sol),sol)
#plotvscos(sol,256,1)
#OVSNC()

#print(sum(np.array(solss)-np.array(solsf)))
