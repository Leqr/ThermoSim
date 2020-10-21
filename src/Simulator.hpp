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
    int pushTimeStep = 10;
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
    
    public:
    std::string getState();
    void simulate();
    
    void solveNonCoupledDiff(bool MMS = false);
    void OVSNonCoupledDiff(double Pe,int n);
    
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
