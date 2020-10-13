//
//  main.cpp
//  CFDSim
//
//  Created by Léonard Equer on 27.09.20.
//

#include <iostream>
#include "Exporter.hpp"
#include "Simulator.hpp"
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>

void testReader(){
    
    std::vector <double> v1 {
        1,
        2,
        3,
        4,
        5
    };
    std::vector <double> v2 {
           1,
           2,
           3,
           45,
           55
       };
    
    std::vector <double> w1 {
        2,
        2,
        3,
        4,
        5
    };
    std::vector <double> w2 {
           2,
           2,
           3,
           45,
           55
    };
    
    Exporter exporter;

    std::vector<std::vector <double>> v{v1,v2};
    std::vector<std::vector <double>> w{w1,w2};

    exporter.fexport(v,w);
    
}

int main(){
    
    double height;
    double diameter;
    int nCells;
    double T0;
    int nCycles;
    int nTimeStepsCycle;
    std::unordered_map<std::string, double> durations;
    double alphaF;
    double alphaS;
    double uf;

    
    //ask the user for the parameters
    std::ifstream setup;
    setup.open("../../src/setup.txt");
    if( !setup.is_open()) {
       std::cerr << "Error: file could not be opened" << std::endl;
       exit(1);
    }
    
    std::string line;
    int k = 1;
    while(getline(setup, line)){
        switch(k){
            case 2:
                height = std::stod(line);
                break;
            case 4:
                diameter = std::stod(line);
                break;
            case 6:
                nCells = std::stoi(line);
                break;
            case 8:
                T0 = std::stod(line);
                break;
            case 10:
                nCycles = std::stoi(line);
                break;
            case 12:
                nTimeStepsCycle = std::stoi(line);
                break;
            case 14:
                durations["charging"] = std::stod(line);
                break;
            case 16:
                durations["idlecd"] = std::stod(line);
                break;
            case 18:
                durations["discharging"] = std::stod(line);
                break;
            case 20:
                durations["idledc"] = std::stod(line);
                break;
            case 22:
                alphaF = std::stod(line);
                break;
            case 24:
                alphaS = std::stod(line);
                break;
            case 26:
                uf = std::stod(line);
                break;
        }
        k += 1;
    }
    
    /*
    std::cout << "Enter the height :" << std::endl;
    std::cin >> height ;
    
    std::cout << "Enter the diameter :" << std::endl;
    std::cin >> diameter ;
    
    std::cout << "Enter the number of cells :" << std::endl;
    std::cin >> nCells ;
    
    std::cout << "Enter the initial temperature :" << std::endl;
    std::cin >> T0 ;
    
    
    std::cout << "Enter the charging state duration :" << std::endl;
    std::cin >> durations["charging"] ;
    
    std::cout << "Enter the idle charging-discharging state duration :" << std::endl;
    std::cin >> durations["idlecd"] ;
    
    std::cout << "Enter the discharging state time :" << std::endl;
    std::cin >> durations["discharging"] ;
    
    std::cout << "Enter the idle discharging-charging state time :" << std::endl;
    std::cin >> durations["idledc"] ;
    
    
    std::cout << "Enter the number of cycles :" << std::endl;
    std::cin >> nCycles ;
    
    std::cout << "Enter the number of time steps per cycle :" << std::endl;
    std::cin >> nTimeStepsCycle ;
    
    std::cout << "Enter alphaF :" << std::endl;
    std::cin >> alphaF ;
    
    std::cout << "Enter alphaS :" << std::endl;
    std::cin >> alphaS ;
    */
    
    //initialize the time serie exporter
    Exporter exporter;

    //run the simulation
    Simulator sim(durations,height,diameter,nCells,T0,nCycles,nTimeStepsCycle,exporter,alphaF,alphaS,uf);
    sim.solveNonCoupledAdvDiffFluid();
    sim.solveNonCoupledDiffSolid();
}
