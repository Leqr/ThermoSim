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

std::string Exporter::getOutFile(){
    return this->outFile;
}

void Exporter::setOutFile(std::string outFile){
    this->outFile = outFile;
}

void Exporter::fexport(std::vector<std::vector<double>> fluidtemp,std::vector<std::vector<double>> solidtemp){
    
    outdata << "fluidtemp" << std::endl;

    for(int i=0; i < int(fluidtemp.size()); i++){
        for(int j=0; j < int(fluidtemp[i].size()); j++){
            outdata << fluidtemp[i][j];
            if(j != int(fluidtemp[i].size())-1){
                outdata << ",";
            }
        }
        outdata << std::endl;
    }
    
    outdata << "solidtemp" << std::endl;
    for(int i=0; i < int(solidtemp.size()); i++){
        for(int j=0; j < int(solidtemp[i].size()); j++){
            outdata << fluidtemp[i][j];
            if(j != int(solidtemp[i].size())-1){
                outdata << ",";
            }
        }
        outdata << std::endl;
    }
}

void Exporter::exportState(std::string state){

    outdatastate << state << std::endl;
}

void Exporter::pushFluid(std::vector<double> fluidtemp){
    for(int j=0; j < int(fluidtemp.size()); j++){
        fluidtempdata << fluidtemp[j];
        if(j != int(fluidtemp.size())-1){
            fluidtempdata << ",";
        }
    }
    fluidtempdata << std::endl;
}

void Exporter::pushSolid(std::vector<double> solidtemp){
    for(int j=0; j < int(solidtemp.size()); j++){
        solidtempdata << solidtemp[j];
        if(j != int(solidtemp.size())-1){
            solidtempdata << ",";
        }
    }
    solidtempdata << std::endl;
}

Exporter::Exporter(std::string outFile,std::string stateFile):outFile(outFile),stateFile(stateFile)
{
    outdata.open(outFile);
    if( !outdata.is_open()) {
       std::cerr << "Error: file could not be opened" << std::endl;
       exit(1);
    }
    
    outdatastate.open(stateFile);
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
}

Exporter::~Exporter()
{
    
    outdata.close();
    outdatastate.close();
    fluidtempdata.close();
    solidtempdata.close();
    
}

