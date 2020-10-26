import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

def plotMov(nCells,nTimeStep,sol):
    #animation of the solution

    fig = plt.figure()
    ax = plt.axes(xlim=(0, nCells),ylim = (200,600))
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

def OVSNC():
    ovs = open("../Data/OVSNonCoupledClean/ovsNCn=1,Pe=0.001.txt")

    ovsdata = []

    for line in ovs:
        a = line.strip().split(',')
        af = []
        for val in a:
            af.append(float(val))
        ovsdata.append(af)
    '''
    fig, axs = plt.subplots(2, 2)

    axs[0, 0].plot([np.log(10/n) for n in [8,16,32,64,128,256,512]],[np.log(val) for val in ovsdata[2]], 'tab:green',marker=".")
    axs[0, 0].set_title('L1 norm for Ts, n = 1')
    axs[0, 0].set_ylabel('Log(E)')
    axs[0, 0].set_xlabel('Log(h)')
    axs[0, 0].set_ylim([-8,5])


    axs[0, 1].plot([np.log(10/n) for n in [8,16,32,64,128,256,512]],[np.log(val) for val in ovsdata[3]], 'tab:orange',marker=".")
    axs[0, 1].set_title('Linf norm for Ts, n = 1')
    axs[0, 1].set_ylabel('Log(E)')
    axs[0, 1].set_xlabel('Log(h)')
    axs[0, 1].set_ylim([-8,5])

    y = [np.log(np.abs((ovsdata[2][i+1]-ovsdata[2][i])/(ovsdata[2][i]-ovsdata[2][i-1])))/np.log(1/2) for i in range(1,len(ovsdata[2])-1)]
    axs[1, 0].plot([np.log(10/n) for n in [16,32,64,128,256]],y ,'tab:green',marker=".")
    axs[1, 0].set_ylabel('p')
    axs[1, 0].set_xlabel('Log(h)')
    axs[1, 0].set_ylim([0,5])

    y2 = [np.log(np.abs((ovsdata[3][i+1]-ovsdata[3][i])/(ovsdata[3][i]-ovsdata[3][i-1])))/np.log(1/2) for i in range(1,len(ovsdata[3])-1)]

    axs[1, 1].plot([np.log(10/n) for n in [16,32,64,128,256]],y2, 'tab:orange',marker=".")
    axs[1, 1].set_ylabel('p')
    axs[1, 1].set_xlabel('Log(h)')
    axs[1, 1].set_ylim([0,5])


    '''
    fig, axs = plt.subplots(2, 2)

    axs[0, 0].plot([np.log(10/n) for n in [8,16,32,64,128,256,512]],[np.log(val) for val in ovsdata[0]], 'tab:green',marker=".")
    axs[0, 0].set_title('L1 norm for Tf, n = 1, Pe = 1')
    axs[0, 0].set_ylabel('Log(E)')
    axs[0, 0].set_xlabel('Log(h)')
    axs[0, 0].set_ylim([-8,20])

    axs[0, 1].plot([np.log(10/n) for n in [8,16,32,64,128,256,512]],[np.log(val) for val in ovsdata[1]], 'tab:orange',marker=".")
    axs[0, 1].set_title('Linf norm for Tf, n = 1, Pe = 1')
    axs[0, 1].set_ylabel('Log(E)')
    axs[0, 1].set_xlabel('Log(h)')
    axs[0, 1].set_ylim([-8,20])


    y = [np.log(np.abs((ovsdata[0][i+1]-ovsdata[0][i])/(ovsdata[0][i]-ovsdata[0][i-1])))/np.log(1/2) for i in range(1,len(ovsdata[0])-1)]
    axs[1, 0].plot([np.log(10/n) for n in [16,32,64,128,256]],y ,'tab:green',marker=".")
    axs[1, 0].set_ylabel('p')
    axs[1, 0].set_xlabel('Log(h)')
    axs[1, 0].set_ylim([0,4])

    y2 = [np.log(np.abs((ovsdata[1][i+1]-ovsdata[1][i])/(ovsdata[1][i]-ovsdata[1][i-1])))/np.log(1/2) for i in range(1,len(ovsdata[1])-1)]

    axs[1, 1].plot([np.log(10/n) for n in [16,32,64,128,256]],y2, 'tab:orange',marker=".")
    axs[1, 1].set_ylabel('p')
    axs[1, 1].set_xlabel('Log(h)')
    axs[1, 1].set_ylim([0,4])

    plt.legend()
    plt.show()

def plotvscos(sol,ncells,n):
    plt.plot(sol[-1])
    plt.plot([np.cos(2*3.14*n/10*x) for x in np.linspace(0,10,ncells)],label = "real cos")
    plt.suptitle("n=128,Ts")
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

#plotMov(16,len(sol),sol)
plotvscos(sol,32,1)
#OVSNC()
