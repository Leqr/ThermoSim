//
//  Exporter.hpp
//  CFDSim
//
//  Created by Léonard Equer on 27.09.20.
//

#ifndef Exporter_hpp
#define Exporter_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>


class Exporter {
    private:
    std::string outFile;
    std::string stateFile;
    std::ofstream outdata;
    std::ofstream outdatastate;
    std::ofstream fluidtempdata;
    std::ofstream solidtempdata;
    
    public:
    std::string getOutFile();
    void setOutFile(std::string outFile);
    void fexport(std::vector<std::vector<double>> fluidtemp,std::vector<std::vector<double>> solidtemp);
    void pushFluid(std::vector<double> fluidtemp);
    void pushSolid(std::vector<double> solidtemp);
    void exportState(std::string state);
    Exporter(std::string outFile = "outCFD.txt",std::string stateFile = "stateCFD.txt");
    ~Exporter();
};

#endif /* Exporter_hpp */
