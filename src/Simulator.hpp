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
    std::unordered_map<std::string, int> durations;
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
    int pushTimeStep = 10;
    int checkSteadyStateTimeStep = 50;
    double errThreshold = 1e-8;
    double dt;
    double dx; //calculated from user param
    double Lbc; //initialized in the constructor
    double Rbc;
    double sim_uf;
    double maximumThermalEnergy;
    
    //parameters for last part of the project (storage design study)
    double ThermalEnergyBefChar;
    double ThermalEnergyBefDischar;
    std::vector<double> capacityFactors;
    std::vector<double> exergyEfficiencies;
    double fluxdout;
    double fluxdin;
    double fluxcout;
    double fluxcin;

    
    //Method of Manufactured solutions parameters
    double n_fluid = 1.0;
    double n_solid = 1.0;
    double k_fluid;
    double k_solid;
    
    //for checking with an analytic solution
    std::vector<double> analyticDataFluid;
    std::vector<double> analyticDataSolid;
    std::vector<double> analyticDataX;
    
    
    public:
    //main simulation function
    void simulate(bool MMS = false,bool coupled = false);
    
    //Solvers and OVS
    void solveDiff(const std::vector<double> &oldsols,const std::vector<double> &oldsolf,std::vector<double> &sols, std::vector<double> &solf, bool MMS = false,bool coupled = false);
    void OVS(double Pe,int n,bool coupled = false);
    
    
    //Utility functions
    std::string getState();
    void checkStabCond();
    double L1Error(const std::vector<double> &numSol,const std::vector<double> &analySol);
    double LinfError(const std::vector<double> &numSol,const std::vector<double> &analySol);
    double computeStoredThermalEnergy(std::vector<double> &sols, std::vector<double> &solf);
    void computeNextExergyFlux(std::vector<double> &oldsolf,std::vector<double> &solf);
    
    //Constructor
    Simulator(std::unordered_map<std::string, int> durations,
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
    
};

#endif /* Simulator_hpp */
