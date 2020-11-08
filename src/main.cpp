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

int main(){
    
    //parameters for the simulation
    double height;
    double diameter;
    int nCells;
    double T0;
    int nCycles;
    int nTimeStepsCycle;
    std::unordered_map<std::string, int> durations;
    double alphaF;
    double alphaS;
    double uf;
    double hvs;
    double hvf;
    
    //open the setup.txt file
    std::ifstream setup;
    
    
    //setup.open("../../src/setup.txt");
    //setup file for part 5 of the project
    setup.open("../../src/setuppr5.txt");
    
    
    if( !setup.is_open()) {
       std::cerr << "Error: setup file could not be opened" << std::endl;
       exit(1);
    }
    
    //read the parameters
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
                durations["charging"] = std::stoi(line);
                break;
            case 16:
                durations["idlecd"] = std::stoi(line);
                break;
            case 18:
                durations["discharging"] = std::stoi(line);
                break;
            case 20:
                durations["idledc"] = std::stoi(line);
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
            case 28:
                hvs = std::stod(line);
                break;
            case 30:
                hvf = std::stod(line);
                break;
        }
        k += 1;
    }

    //initialize the time serie exporter
    Exporter exporter;

    //run the simulation
    Simulator sim(durations,height,diameter,nCells,T0,nCycles,nTimeStepsCycle,exporter,alphaF,alphaS,uf,hvs,hvf);
    //sim.OVSNonCoupledDiff(1,2);
    sim.simulate(false,true);
   
}

/*
 ******************************************************************************
 Deprecated, for reference
 ******************************************************************************
*/

/*
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
 */
