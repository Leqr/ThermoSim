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
#include <stdlib.h>

int main(int argc, char** argv){

    std::string pathToSrc("../src/");

    //mode setup from executable call
    int mode = std::stoi(argv[1]);

    //open the setup.txt file
    std::ifstream setup;

    switch(mode){
        case 0:
            //setup file for MMS
            setup.open("../src/setupMMS.txt");
            break;
        case 1:
            //setup file for non coupled OVS
            setup.open(pathToSrc + "setupOVS.txt");
            break;
        case 2:
            //setup file for coupled OVS
            setup.open(pathToSrc + "setupOVS.txt");
            break;
        case 3:
            //setup file for part 5 of the project
            setup.open(pathToSrc + "setuppr5.txt");
            break;
        case 4:
            //setup for last part, total simulation
            setup.open(pathToSrc + "setup.txt");
            break;
    }


    if(!setup.is_open()) {
      std::cerr << "Error: setup file could not be opened" << std::endl;
      exit(1);
    }

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
    Exporter exporter(pathToSrc);

    //initialize the simulator
    Simulator sim(durations,height,diameter,nCells,T0,nCycles,nTimeStepsCycle,exporter,alphaF,alphaS,uf,hvs,hvf);
    std::string cmd;
    switch(mode){
        case 0:
            //solve the uncoupled equation with MMS source term and plot it
            std::cout << "Numerical solution plotting against manufactured solution." << std::endl;
            sim.simulate(true,false);
            cmd = "python " + pathToSrc + "plotter.py 0";
            system(cmd.c_str());
            break;

        case 1:
            //non coupled OVS and plot
            std::cout << "OVS Non Coupled." << std::endl;
            sim.OVS(1000,1,false);
            cmd = "python " + pathToSrc + "plotter.py 1";
            system(cmd.c_str());
            break;

        case 2:
            //coupled OVS and plot
            std::cout << "OVS Coupled." << std::endl;
            sim.OVS(1000,1,true);
            cmd = "python " + pathToSrc + "plotter.py 2";
            system(cmd.c_str());
            break;

        case 3:
            //solve the coupled equation for part 5 and plots the result
            std::cout << "Analytic vs numeric." << std::endl;
            sim.simulate(false,true);
            cmd = "python " + pathToSrc + "plotter.py 3";
            system(cmd.c_str());
            break;

        case 4:
            //solve the coupled equation for design study and plots the result + show the exergy and capacity in the terminal
            std::cout << "Design Study" << std::endl;
            std::cout << "Diameter : " << diameter << std::endl;
            sim.simulate(false,true);
            cmd = "python " + pathToSrc + "plotter.py 4";
            system(cmd.c_str());
            break;

    }
    return 0;
}
