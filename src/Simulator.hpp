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
    
    //sim parameters
    int pushTimeStep = 100;
    int checkSteadyStateTimeStep = 1000;
    double errThreshold = 1e-6;
    double dt = 0.5;
    double dx;
    double Lbc;
    
    //Method of Manufactured solutions
    double n = 1.0;
    double k = 2*M_PI*n/height;
    
    
    public:
    std::string getState();
    void simulate(bool MMS = false);
    
    void solveNonCoupledDiff(bool MMS = false);
    std::vector<std::vector<double>> solveDiff(std::vector<double> oldsols, std::vector<double> oldsolf, bool MMS = false,bool coupled = true);
    void OVSNonCoupledDiff(double Pe,int n);
    
    void checkStabCond();
    double L1Error(std::vector<double> numSol,std::vector<double> analySol);
    double LinfError(std::vector<double> numSol,std::vector<double> analySol);
    
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
              double uf);
    
};

#endif /* Simulator_hpp */
