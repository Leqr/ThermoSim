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

    std::ofstream outdatastate;
    std::ofstream fluidtempdata;
    std::ofstream solidtempdata;
    std::ofstream ovsdata;
    std::string pathToSrc;

    public:

    void pushFluid(const std::vector<double> &fluidtemp);
    void pushSolid(const std::vector<double> &solidtemp);
    void pushOVS(const std::vector<double> &l);
    void exportState(std::string state);
    Exporter(std::string pathToSrc);
    ~Exporter();

};

#endif /* Exporter_hpp */
