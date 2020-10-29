//
//  Simulator.hpp
//  CFDSim
//
//  Created by Léonard Equer on 02.10.20.
//

#ifndef Simulator_hpp
#define Simulator_hpp

#include <stdio.h>
#include <string>
#include <unordered_map>
#include "Exporter.hpp"


class Simulator {
    private:
    
    // initialised by the user
    std::unordered_map<std::string, double> durations;
    double height;
    double diameter;
    int nCells;
    double T0;
    int nCycles;
    std::string state;
    Exporter& exporter;
    int nTimeStepsCycle;
    double alphaF;
    double alphaS;
    double uf;
    double hvf;
    double hvs;
    
    //simulation parameters
    int pushTimeStep = 1000;
    int checkSteadyStateTimeStep = 50;
    double errThreshold = 1e-6;
    double dt = 1;
    double dx; //calculated from user param
    double Lbc; //initialized in the constructor
    
    //Method of Manufactured solutions parameters
    double n_fluid = 1.0;
    double n_solid = 1.0;
    double k_fluid;
    double k_solid;
    
    
    public:
    //main simulation function
    void simulate(bool MMS = false,bool coupled = false);
    
    //Solvers and OVS
    void solveDiff(const std::vector<double> &oldsols,const std::vector<double> &oldsolf,std::vector<double> &sols, std::vector<double> &solf, bool MMS = false,bool coupled = false);
    void OVSNonCoupledDiff(double Pe,int n);
    
    //Utility functions
    std::string getState();
    void checkStabCond();
    double L1Error(const std::vector<double> &numSol,const std::vector<double> &analySol);
    double LinfError(const std::vector<double> &numSol,const std::vector<double> &analySol);
    
    //Constructor
    Simulator(std::unordered_map<std::string, double> durations,
              double height,
              double diameter,
              int nCells,
              double T0,
              int nCycles,
              int nTimeStepsCycle,
              Exporter& exporter,
              double alphaF,
              double alphaS,
              double uf,
              double hvs,
              double hvf);
    
    /******************** Deprecated ******************/
    void solveNonCoupledDiff(bool MMS = false);

    
};

#endif /* Simulator_hpp */
