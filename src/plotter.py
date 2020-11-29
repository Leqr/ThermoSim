#
#  Plotter.py
#  CFDSim
#
#  Created by Léonard Equer on 02.10.20.
#

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from matplotlib.animation import FuncAnimation
import sys
import os

def importRealSols(pathToSrc):
    #loads the analytical solution
    s = open(pathToSrc+"sol-exact-5.00000E+03.dat")

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

def importComputedSols(pathToBuild):
    #loading the solutions
    ffluid = open(pathToBuild+"fluid.txt")

    solsf = []

    for line in ffluid:
        a = line.strip().split(',')
        af = []
        for val in a:
            af.append(float(val))
        solsf.append(af)

    fsolid = open(pathToBuild+"solid.txt")

    solss = []

    for line in fsolid:
        a = line.strip().split(',')
        af = []
        for val in a:
            af.append(float(val))
        solss.append(af)

    return (solsf,solss)


def readState(pathToBuild):
    #read the state trace of the program
    s = open(pathToBuild+"stateCFD.txt")

    statedata = s.readlines()

    plt.plot(statedata)
    plt.show()


def plotMov(nCells,nTimeStep,sol):
    #animation of one solution
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
                                   frames=nTimeStep, interval=5,blit = False)
    plt.show()

def plotMovDouble(nCells,nTimeStep,sol):
    #animation of the solution with both fluid and solid
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
                                   frames=nTimeStep, interval=5, blit = False)
    plt.legend()
    plt.show()

def OVSNC(pathToBuild):
    #plot the ovs results
    ovs = open(pathToBuild+"ovs.txt")
    ovsdata = []
    for line in ovs:
        a = line.strip().split(',')
        af = []
        for val in a:
            af.append(float(val))
        ovsdata.append(af)

    #generate the plots
    fig, axs = plt.subplots(2, 2)
    vals = [16,32,64,128,256]

    axs[0, 0].plot([np.log(1.0/n) for n in vals],[np.log(val) for val in ovsdata[2]],marker=".",label = "L1 error")
    axs[0, 0].set_title('Tsolid')
    axs[0, 0].set_ylabel('Log(E)')
    axs[0, 0].set_xlabel('Log(h)')
    axs[0, 0].set_ylim([-8,10])
    axs[0, 0].plot([np.log(1.0/n) for n in vals],[np.log(val) for val in ovsdata[3]],marker=".",label = "Linf error")
    axs[0, 0].legend(loc="upper right")

    y = [np.log(ovsdata[2][i]/ovsdata[2][i-1])/np.log(0.5) for i in range(1,len(ovsdata[2]))]
    axs[1, 0].plot([np.log(1.0/n) for n in vals[1:]],y ,marker=".",label = "L1")
    axs[1, 0].set_ylabel('p')
    axs[1, 0].set_xlabel('Log(h)')
    axs[1, 0].set_ylim([0,4])
    y2 = [np.log(ovsdata[3][i]/ovsdata[3][i-1])/np.log(0.5) for i in range(1,len(ovsdata[3]))]
    axs[1, 0].plot([np.log(1.0/n) for n in vals[1:]],y2,marker=".",label = "Linf")
    axs[1, 0].legend(loc="upper right")


    axs[0, 1].plot([np.log(1.0/n) for n in vals],[np.log(val) for val in ovsdata[0]],marker=".",label = "L1 error")
    axs[0, 1].set_title('Tfluid')
    axs[0, 1].set_ylabel('Log(E)')
    axs[0, 1].set_xlabel('Log(h)')
    axs[0, 1].set_ylim([-8,10])
    axs[0, 1].plot([np.log(1.0/n) for n in vals],[np.log(val) for val in ovsdata[1]],marker=".",label = "Linf error")
    axs[0, 1].legend(loc="upper right")

    y = [np.log(ovsdata[0][i]/ovsdata[0][i-1])/np.log(0.5) for i in range(1,len(ovsdata[0]))]
    axs[1, 1].plot([np.log(1.0/n) for n in vals[1:]],y ,marker=".",label = "L1")
    axs[1, 1].set_ylabel('p')
    axs[1, 1].set_xlabel('Log(h)')
    axs[1, 1].set_ylim([0,4])
    y2 = [np.log(ovsdata[1][i]/ovsdata[1][i-1])/np.log(0.5) for i in range(1,len(ovsdata[1]))]
    axs[1, 1].plot([np.log(1.0/n) for n in vals[1:]],y2,marker=".",label = "Linf")
    axs[1, 1].legend(loc="upper right")

    plt.legend()
    plt.show()

def plotvscos(sol,ncells,n):
    #plot against the MMS cosine
    plt.plot(sol[-1],label = "numerical solution")
    plt.plot([np.cos(2*3.14*n/10*x) for x in np.linspace(0,10,ncells)],label = "real cos")
    plt.suptitle("Ts,nC=256,n=1")
    plt.legend()
    plt.show()

def plotvsreal(computed,real):
    #plot against the analytical solution
    sol_fluid,sol_solid = computed
    x,real_sol_fluid,real_sol_solid = real

    plt.plot(x,sol_fluid[-1],label = 'Numerical Tf')
    plt.plot(x,real_sol_fluid, label = 'Analytic Tf')
    plt.plot(x,sol_solid[-1],label = 'Numerical Ts')
    plt.plot(x,real_sol_solid,label = 'Analytic Ts')
    plt.ylabel('Celsius')
    plt.xlabel('Length [m]')

    plt.title("Numerical solution against analytical solution.")
    plt.legend()
    plt.show()

def getH_fFromD(d):
    #parameters calculation for each d for design study
    uf = 10/(np.pi*1835*np.power(d/2,2)*0.4)
    Re = 0.4*1835.6*uf*0.03/2.63
    Pr = 2.63*1511/0.52
    Nufs = (0.255/0.4)*np.power(Pr,1/3)*np.power(Re,2/3)
    hfs = Nufs*(0.52/0.03)
    h = 1/((1/hfs)+(0.03/(10*2.0)))
    hv = h*6*0.6/0.03
    hvf = hv/(0.4*1835.6*1511.8)
    hvs = hv/(0.6*2600*900)
    height = 300/(np.pi*np.power(d/2,2))
    print("uf : " + str(uf) + ", height : " + str(height) + ", hvf : " + str(hvf) + ", hvs : "+str(hvs))



if __name__ == "__main__":

    pathToBuild = ""
    pathToSrc = ""

    mode = -1
    try :
        mode = int(sys.argv[1])
        pathToBuild = ""
        pathToSrc = "../src/"

    except :
        # no argument passed, means manual run of the script for debugging
        mode = -1
        pathToBuild = "../build/"
        pathToSrc = ""

    sol_fluid,sol_solid = importComputedSols(pathToBuild)

    sol = sol_fluid

    #used for debugging and running directly python plotter.py without argv
    '''Plots a Movie with the two computed solutions'''
    #nCells = 128
    #plotMovDouble(nCells,len(sol_solid),[sol_fluid,sol_solid])

    '''Plots a Movie with one computed solution'''
    #plotMov(nCells,len(sol),sol)

    '''Plot the last computed solution against a cosine for MMS '''
    #n = 4
    #nCells = 64
    #plotvscos(sol,nCells,n)

    '''Plots the OVS results'''
    #OVSNC(pathToBuild)

    '''Plots the state trace of the simulation'''
    #readState(pathToBuild)

    '''Plots the analytical solution against the computed sol'''
    #x,real_sol_fluid,real_sol_solid = importRealSols(pathToSrc)
    #nCells = 1024
    #plotvsreal((sol_fluid,sol_solid),(x,real_sol_fluid,real_sol_solid))

    '''Calculate the hvs and hvf values depending on D'''
    '''
    for d in [4.0,5.0,6.0,7.0,8.0]:
        getH_fFromD(d)
    '''

    #mode ran from c++ code
    if mode == 0 :
        '''Plot the last computed solution against a cosine for MMS '''
        n = 1
        nCells = 64
        plotvscos(sol,nCells,n)
    if mode == 1 :
        '''Plots the OVS results'''
        OVSNC(pathToBuild)
    if mode == 2 :
        '''Plots the OVS results'''
        OVSNC(pathToBuild)
    if mode == 3:
        '''Plots the analytical solution against the computed sol'''
        x,real_sol_fluid,real_sol_solid = importRealSols(pathToSrc)
        nCells = 1024
        plotvsreal((sol_fluid,sol_solid),(x,real_sol_fluid,real_sol_solid))
    if mode == 4:
        '''Plots a Movie with the two computed solutions'''
        nCells = 128
        plotMovDouble(nCells,len(sol_solid),[sol_fluid,sol_solid])
