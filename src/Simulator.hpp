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
    std::unordered_map<std::string, double> durations;
    double height;
    double diameter;
    int nCells;
    double T0;
    int nCycles;
    int nTimeStepsCycle;
    std::string state;
    Exporter& exporter;
    double alphaF;
    double alphaS;
    double uf;
    
    public:
    std::string getState();
    void simulate();
    void solveNonCoupledAdvDiffFluid();
    void solveNonCoupledDiffSolid();
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
