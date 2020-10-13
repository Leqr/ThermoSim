//
//  Simulator.cpp
//  CFDSim
//
//  Created by Léonard Equer on 02.10.20.
//

#include "Simulator.hpp"
#include <unordered_map>
#include <iostream>

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

void Simulator::solveNonCoupledAdvDiffFluid(){
    double T = 100.0;
    double dt = T/nTimeStepsCycle;
    double dx = 1/double(nCells);
    //initial temperature
    std::vector<double> sol(nCells,T0);
    
    //timestep simulation of the program
    for(int i = 1; i <= nTimeStepsCycle; ++i){
        
        exporter.pushFluid(sol);
        
        std::vector<double> oldsol = sol;
        sol.clear();
        
        //over all space for each timestep
        for(int j = 0; j<=nCells-1;++j){
            if (j==0){
                double val;
                val = 340 - (uf*dt/dx)*(oldsol[j]) + (alphaF*dt/(dx*dx))*(oldsol[j+1]-2*oldsol[j]);
                sol.push_back(val);
            }
            if (j > 0 and j < (nCells-1)){
                double val;
                val = oldsol[j] - (uf*dt/dx)*(oldsol[j]-oldsol[j-1]) + (alphaF*dt/(dx*dx))*(oldsol[j+1]-2*oldsol[j]+ oldsol[j-1]);
                sol.push_back(val);
            }
            if (j==nCells-1){
                double val;
                val = oldsol[j] - (uf*dt/dx)*(oldsol[j]-oldsol[j-1]) + (alphaF*dt/(dx*dx))*(oldsol[j-1]-2*oldsol[j]);
                sol.push_back(val);
            }
        }
    }
}

void Simulator::solveNonCoupledDiffSolid(){
    double T = 1000.0;
    double dt = T/nTimeStepsCycle;
    double dx = 1/double(nCells);
    //initial temperature
    std::vector<double> sol(nCells,T0);
    
    //timestep simulation of the program
    for(int i = 1; i <= nTimeStepsCycle; ++i){
        
        exporter.pushSolid(sol);
        
        std::vector<double> oldsol = sol;
        sol.clear();
        
        //over all space for each timestep
        for(int j = 0; j<=nCells-1;++j){
            if (j==0){
                double val;
                //val = oldsol[j] + (alphaS*dt/(dx*dx))*(oldsol[j+1]-2*oldsol[j]);
                val = 400 + (alphaS*dt/(dx*dx))*(oldsol[j+1]-2*oldsol[j]);

                sol.push_back(val);
            }
            if (j > 0 and j < (nCells-1)){
                double val;
                val = oldsol[j] + (alphaS*dt/(dx*dx))*(oldsol[j+1]-2*oldsol[j]+ oldsol[j-1]);
                sol.push_back(val);
            }
            if (j==nCells-1){
                double val;
                val = 293 + (alphaF*dt/(dx*dx))*(oldsol[j-1]-2*293);
                sol.push_back(val);
            }
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
