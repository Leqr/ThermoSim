import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from matplotlib.animation import FuncAnimation

def importRealSols():
    s = open("sol-exact-5.00000E+03.dat")

    x = []
    fluid = []
    solid = []

    i = 0
    for line in s:
        if i != 0:
            a = line.strip().split(' ')
            af = []
            for val in a:
                if val != '':
                    af.append(float(val))
            x.append(af[0])
            fluid.append(af[1])
            solid.append(af[2])
        i += 1

    return (x,fluid,solid)

def importComputedSols():
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

    return (solsf,solss)


def readState():
    s = open("../build/debug/stateCFD.txt")

    statedata = s.readlines()

    plt.plot(statedata)
    plt.show()


def plotMov(nCells,nTimeStep,sol):
    #animation of the solution

    fig = plt.figure()
    ax = plt.axes(xlim=(0, nCells-1),ylim = (200,1000))
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
    ln0, = plt.plot([], [], 'cornflowerblue',label = "Fluid Temperature")
    ln1, = plt.plot([], [], 'tan',label = "Solid Temperature")

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

def OVSNC(plotsolid):
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
    if plotsolid:
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


    else:
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

def plotvsreal(computed,real):
    sol_fluid,sol_solid = computed
    x,real_sol_fluid,real_sol_solid = real

    plt.plot(x,sol_fluid[-1],label = 'Numerical Tf')
    plt.plot(x,real_sol_fluid, label = 'Analytic Tf')
    plt.plot(x,sol_solid[-1],label = 'Numerical Ts')
    plt.plot(x,real_sol_solid,label = 'Analytic Ts')

    plt.title("Numerical solution against analytical solution.")
    plt.legend()
    plt.show()

if __name__ == "__main__":

    nCells = 1024

    sol_fluid,sol_solid = importComputedSols()
    x,real_sol_fluid,real_sol_solid = importRealSols()

    '''Plots a Movie with the two computed solutions'''
    #plotMovDouble(nCells,len(sol_solid),[sol_fluid,sol_solid])

    sol = sol_fluid

    '''Plots a Movie with one computed solution'''
    #plotMov(nCells,len(sol),sol)

    '''Plot the last computed solution against a cosine for MMS '''
    k = 1
    #plotvscos(sol,nCells,k)

    '''Plots the OVS results'''
    #OVSNC(false)

    '''Plots the state trace of the simulation'''
    #readState()

    '''Plots the analytical solution against the computed sol'''
    nCells = 1024
    plotvsreal((sol_fluid,sol_solid),(x,real_sol_fluid,real_sol_solid))
