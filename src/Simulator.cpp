//
//  Simulator.cpp
//  CFDSim
//
//  Created by Léonard Equer on 02.10.20.
//

#include "Simulator.hpp"
#include <unordered_map>
#include <iostream>
#import <math.h>

std::string Simulator::getState(){
    return state;
}

void Simulator::simulate(){
    
    //timestep simulation of the program
    for(int i = 1; i < nTimeStepsCycle*nCycles; ++i){
        int val = i%nTimeStepsCycle;
        //updating the state of the simulation (might find a way more efficient but less clean)
        if (val == durations["charging"]){
            state = "charging";
        }
        if (val == durations["idlecd"]+durations["charging"]){
            state = "idlecd";
        }
        if (val == durations["discharging"]+durations["idlecd"]+durations["charging"]){
            state = "discharging";
        }
        if (val == 0){
            state = "idledc";
        }
        exporter.exportState(state);
        
    }
}


void Simulator::solveNonCoupledDiff(bool MMS){
    
    //param
    double T = 5000.0;
    double dt = T/nTimeStepsCycle;
    
    double dx = height/double(nCells);
    
    double Lbc = 340;
    double Rbc = 293;
    
    double k = 2*M_PI*1/height;
    
    double errThreshold = 1e-3;
    
    //initial temperature
    std::vector<double> solf(nCells,T0);
    std::vector<double> sols(nCells,T0);
    
    //to have these variables as global
    std::vector<double> oldsols = sols;
    std::vector<double> oldsolf = solf;
    
    //source term for Method of Manufactured solutions
    double source = 0;
    
    double sigma = uf*dt/dx;
    double df = alphaF*dt/(dx*dx);
    double ds = alphaS*dt/(dx*dx);
    
    if (2*ds > 1){
        std::cerr << "ds is too big, the method might not be stable"  << std::endl;
    }
    
    if((sigma*sigma > sigma + 2*df) and (sigma + 2*df > 1)){
        std::cerr << "The method might not be stable"  << std::endl;
    }
    
    
    
    //timestep simulation of the program
    for(int i = 1; i <= nTimeStepsCycle; ++i){
        
        if (i%this->pushTimeStep == 0){
            exporter.pushSolid(sols);
            exporter.pushFluid(solf);
        }
        
        if(MMS){
            source = alphaS*k*k*cos(k*dx*i);
        }
        
        std::vector<double> oldsols = sols;
        std::vector<double> oldsolf = solf;

        sols.clear();
        solf.clear();
        
        //over all space for each timestep
        
        //solid
        for(int j = 0; j<=nCells-1;++j){
            if (j==0){
                double val;
                //val = oldsol[j] + (alphaS*dt/(dx*dx))*(oldsol[j+1]-2*oldsol[j]);
                val = Lbc + (alphaS*dt/(dx*dx))*(oldsols[j+1]-2*oldsols[j])+dt*source;

                sols.push_back(val);
            }
            if (j > 0 and j < (nCells-1)){
                double val;
                val = oldsols[j] + (alphaS*dt/(dx*dx))*(oldsols[j+1]-2*oldsols[j]+ oldsols[j-1]) + dt*source;
                sols.push_back(val);
            }
            if (j==nCells-1){
                double val;
                val = Rbc + (alphaF*dt/(dx*dx))*(oldsols[j-1]-2*Rbc) + dt*source;
                sols.push_back(val);
            }
            
            //fluid
            if (j==0){
                double val;
                val = Lbc - (uf*dt/dx)*(oldsolf[j]) + (alphaF*dt/(dx*dx))*(oldsolf[j+1]-2*oldsolf[j])+dt*source;
                solf.push_back(val);
            }
            if (j > 0 and j < (nCells-1)){
                double val;
                val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-oldsolf[j-1]) + (alphaF*dt/(dx*dx))*(oldsolf[j+1]-2*oldsolf[j]+ oldsolf[j-1])+dt*source;
                solf.push_back(val);
            }
            if (j==nCells-1){
                double val;
                val = Rbc - (uf*dt/dx)*(oldsolf[j]-oldsolf[j-1]) + (alphaF*dt/(dx*dx))*(oldsolf[j-1]-2*Rbc)+dt*source;
                solf.push_back(val);
            }
        }
    }
    if(MMS){
        int err = 0;
        for(int i = 0; i < int(solf.size());++i){
            err += 1/dt * abs(solf[i]-oldsolf[i]);
        }
        err = 1.0/nCells*err;
        if(err > errThreshold){
            std::cerr << "Steady state not attained for threshold " << errThreshold << std::endl;
        }
    }
}

//starts the simulation in the charging state
Simulator::Simulator(std::unordered_map<std::string, double> durations,
                     double height,
                     double diameter,
                     int nCells,
                     double T0,
                     int nCycles,
                     int nTimeStepsCycle,
                     Exporter& exporter,
                     double alphaF,
                     double alphaS,
                     double uf):durations(durations),height(height),diameter(diameter),
                    nCells(nCells),T0(T0),nCycles(nCycles),exporter(exporter),nTimeStepsCycle(nTimeStepsCycle),alphaF(alphaF),alphaS(alphaS),uf(uf),state("idledc")
{}
