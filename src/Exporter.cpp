//
//  Exporter.cpp
//  CFDSim
//
//  Created by Léonard Equer on 27.09.20.
//

#include "Exporter.hpp"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>


void Exporter::exportState(std::string state){

    outdatastate << state << std::endl;
}

void Exporter::pushFluid(const std::vector<double> &fluidtemp){
    for(int j=0; j < int(fluidtemp.size()); j++){
        fluidtempdata << fluidtemp[j];
        if(j != int(fluidtemp.size())-1){
            fluidtempdata << ",";
        }
    }
    fluidtempdata << std::endl;
}

void Exporter::pushSolid(const std::vector<double> &solidtemp){
    for(int j=0; j < int(solidtemp.size()); j++){
        solidtempdata << solidtemp[j];
        if(j != int(solidtemp.size())-1){
            solidtempdata << ",";
        }
    }
    solidtempdata << std::endl;
}

void Exporter::pushOVS(const std::vector<double> &l){
    for(int j=0; j < int(l.size()); j++){
        ovsdata << l[j];
        if(j != int(l.size())-1){
            ovsdata << ",";
        }
    }
    ovsdata << std::endl;
}

Exporter::Exporter(std::string pathToSrc):pathToSrc(pathToSrc)
{

    outdatastate.open("stateCFD.txt");
    if( !outdatastate.is_open()) {
       std::cerr << "Error: file could not be opened" << std::endl;
       exit(1);
    }

    fluidtempdata.open("fluid.txt");
    if( !fluidtempdata.is_open()) {
       std::cerr << "Error: file could not be opened" << std::endl;
       exit(1);
    }

    solidtempdata.open("solid.txt");
    if( !solidtempdata.is_open()) {
       std::cerr << "Error: file could not be opened" << std::endl;
       exit(1);
    }

    ovsdata.open("ovs.txt");
    if( !ovsdata.is_open()) {
       std::cerr << "Error: file could not be opened" << std::endl;
       exit(1);
    }
}

Exporter::~Exporter()
{

    outdatastate.close();
    fluidtempdata.close();
    solidtempdata.close();
    ovsdata.close();

}
