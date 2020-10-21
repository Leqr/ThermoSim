//
//  Simulator.cpp
//  CFDSim
//
//  Created by Léonard Equer on 02.10.20.
//

#include "Simulator.hpp"
#include <unordered_map>
#include <iostream>
#import <math.h>

std::string Simulator::getState(){
    return state;
}

void Simulator::simulate(){
    
    //timestep simulation of the program
    for(int i = 1; i < nTimeStepsCycle*nCycles; ++i){
        int val = i%nTimeStepsCycle;
        //updating the state of the simulation (might find a way more efficient but less clean)
        if (val == durations["charging"]){
            state = "charging";
        }
        if (val == durations["idlecd"]+durations["charging"]){
            state = "idlecd";
        }
        if (val == durations["discharging"]+durations["idlecd"]+durations["charging"]){
            state = "discharging";
        }
        if (val == 0){
            state = "idledc";
        }
        exporter.exportState(state);
        
    }
}


void Simulator::solveNonCoupledDiff(bool MMS){
    
    //param
    double dt = 0.05;
    double dx = height/double(nCells);
    
    std::cout << "dt : " << dt << std::endl;
    std::cout << "dx : " << dx << std::endl;

    //boundary conditions
    double Lbc = 340;
    double Rbc = 293;
        
    //error threshold defining the steady state solution
    double errThreshold = 1e-6;
    //bool of ss attained, when both true the program stops
    bool ss_fluid = false;
    bool ss_solid = false;
    
    //initial temperature
    std::vector<double> solf(nCells,T0);
    std::vector<double> sols(nCells,T0);
    
    //to have these variables as global
    std::vector<double> oldsols = sols;
    std::vector<double> oldsolf = solf;
    
    
    //Method of Manufactured solutions
    double n = 1.0;
    double k = 2*M_PI*n/height;
    double sources = 0;
    double sourcef = 0;
    //puts the right boundary conditions for the MMS
    if(MMS){
        Lbc = 1;
        Rbc = 1;
    }
    
    //check stability conditions
    
    double sigma = uf*dt/dx;
    double df = alphaF*dt/(dx*dx);
    double ds = alphaS*dt/(dx*dx);
    std::cout << "df : " << df << std::endl;
    std::cout << "ds : " << ds << std::endl;
    std::cout << "sigma : " << sigma << std::endl;
    std::cout << std::endl;
    
    if (2*ds > 1){
        std::cerr << "ds is too big, the method might not be stable"  << std::endl;
    }
    
    if(sigma*sigma > sigma + 2*df){
        std::cerr << "The method might not be stable (lower bound)"  << std::endl;
        std::cout <<sigma*sigma << " " <<sigma + 2*df << std::endl;
    }
    
    if(sigma + 2*df > 1){
        std::cerr << "The method might not be stable (upper bound)"  << std::endl;
    }
    
    
    
    //timestep simulation of the program
    for(int i = 1; i <= nTimeStepsCycle; ++i){
        
        //export every multiple of pushtimestep solution to a file
        if (i%this->pushTimeStep == 0){
            exporter.pushSolid(sols);
            exporter.pushFluid(solf);
        }
        
        std::vector<double> oldsols = sols;
        std::vector<double> oldsolf = solf;

        sols.clear();
        solf.clear();
        
        //over all space for each timestep

        //solid
        for(int j = 0; j<=nCells-1;++j){
            
            //calculate the source term in the MMS method
            if(MMS){
                double xi = dx*j;
                sources = alphaS*k*k*cos(k*xi);
                sourcef = alphaF*k*k*cos(k*xi) - uf*k*sin(k*xi);
            }
            
            //solid
            if (j==0){
                double val;
                //val = oldsol[j] + (alphaS*dt/(dx*dx))*(oldsol[j+1]-2*oldsol[j]);
                val = oldsols[j] + (alphaS*dt/(dx*dx))*(oldsols[j+1]-2*oldsols[j] + Lbc)+dt*sources;

                sols.push_back(val);
            }
            if (j > 0 and j < (nCells-1)){
                double val;
                val = oldsols[j] + (alphaS*dt/(dx*dx))*(oldsols[j+1]-2*oldsols[j]+ oldsols[j-1]) + dt*sources;
                sols.push_back(val);
            }
            if (j==nCells-1){
                double val;
                val = oldsols[j] + (alphaF*dt/(dx*dx))*(oldsols[j-1]-2*oldsols[j]+ Rbc) + dt*sources;
                sols.push_back(val);
            }
            
            //fluid
            if (j==0){
                double val;
                val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-Lbc) + (alphaF*dt/(dx*dx))*(oldsolf[j+1]-2*oldsolf[j]+Lbc)+dt*sourcef;
                solf.push_back(val);
            }
            if (j > 0 and j < (nCells-1)){
                double val;
                val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-oldsolf[j-1]) + (alphaF*dt/(dx*dx))*(oldsolf[j+1]-2*oldsolf[j]+ oldsolf[j-1])+dt*sourcef;
                solf.push_back(val);
            }
            if (j==nCells-1){
                double val;
                val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-oldsolf[j-1]) + (alphaF*dt/(dx*dx))*(oldsolf[j-1]-2*oldsolf[j]+Rbc)+dt*sourcef;
                solf.push_back(val);
            }
        }
        
        //check if steady state is attained
        if(MMS){
            
            if (ss_fluid == false){
                double err = 0;
                for(int i = 0; i <nCells;++i){
                    err += 1/dt * abs(solf[i]-oldsolf[i]);
                }
                err = 1.0/nCells*err;
                if(err < errThreshold){
                    std::cerr << "Steady state of fluid solution attained for threshold " << errThreshold << " after " << i << " timesteps." << std::endl;
                    ss_fluid = true;
                }
            }
            
            if (ss_solid == false){
                double err = 0;
                for(int i = 0; i < nCells;++i){
                    err += 1/dt * abs(sols[i]-oldsols[i]);
                }
                err = 1.0/nCells*err;
                if(err < errThreshold){
                    std::cerr << "Steady state of solid solution attained for threshold " << errThreshold << " after " << i << " timesteps." << std::endl;
                    ss_solid = true;
                }
            }
        }
        
        //break out of the time step loop when steady state is attained for both solutions
        if(ss_fluid and ss_solid){
            break;
        }
    }
    if(MMS){
        
        std::vector<double> analyticSol;
        for(int i = 0; i < nCells; ++i){
            analyticSol.push_back(cos(k*dx*i));
        }

        std::cout <<"L1 error : " << L1Error(solf,analyticSol) << std::endl;
        std::cout << "Linf error : " << LinfError(solf,analyticSol) << std::endl;
    }
}

void Simulator::OVSNonCoupledDiff(double Pe, int n){
    double dt = 0.1;
    double alphaf = alphaF;
    
    std::vector<double> L1F;
    std::vector<double> LinfF;
    std::vector<double> L1S;
    std::vector<double> LinfS;
    
    int NCellsi = 8;
    int NCells = 8;
    
    for(int i = 0;i < 7; ++i){
        
        NCells = NCellsi*pow(2,i);
        
        std::cout << "Number of cells : " << NCells << std::endl;
        //param
        double dx = height/double(NCells);
        alphaf= (uf*dx*NCells)/Pe;

        std::cout << "dt : " << dt << std::endl;
        std::cout << "dx : " << dx << std::endl;

        
        //error threshold defining the steady state solution
        double errThreshold = 1e-6;
        //bool of ss attained, when both true the program stops
        bool ss_fluid = false;
        bool ss_solid = false;
        
        //initial temperature
        std::vector<double> solf(NCells,T0);
        std::vector<double> sols(NCells,T0);
        
        //to have these variables as global
        std::vector<double> oldsols = sols;
        std::vector<double> oldsolf = solf;
        
        
        //Method of Manufactured solutions
        double k = 2*M_PI*n/height;
        double sources = 0;
        double sourcef = 0;
        //puts the right boundary conditions for the MMS
        double Lbc = 1;
        double Rbc = 1;
        
        //check stability conditions
        
        double sigma = uf*dt/dx;
        double df = alphaf*dt/(dx*dx);
        double ds = alphaS*dt/(dx*dx);
        std::cout << "df : " << df << std::endl;
        std::cout << "ds : " << ds << std::endl;
        std::cout << "sigma : " << sigma << std::endl;
        std::cout << std::endl;
        
        if (2*ds > 1){
            std::cerr << "ds is too big, the method might not be stable"  << std::endl;
        }
        
        if(sigma*sigma > sigma + 2*df){
            std::cerr << "The method might not be stable (lower bound)"  << std::endl;
            std::cout <<sigma*sigma << " " <<sigma + 2*df << std::endl;
        }
        
        if(sigma + 2*df > 1){
            std::cerr << "The method might not be stable (upper bound)"  << std::endl;
        }
        
        
        
        //timestep simulation of the program
        for(int i = 1; i <= nTimeStepsCycle; ++i){
            
            
            std::vector<double> oldsols = sols;
            std::vector<double> oldsolf = solf;

            sols.clear();
            solf.clear();
            
            //over all space for each timestep

            //solid
            for(int j = 0; j<NCells;++j){
                
                //calculate the source term in the MMS method
                double xi = dx*j;
                sources = alphaS*k*k*cos(k*xi);
                sourcef = alphaf*k*k*cos(k*xi) - uf*k*sin(k*xi);
                
                //solid
                if (j==0){
                    double val;
                    //val = oldsol[j] + (alphaS*dt/(dx*dx))*(oldsol[j+1]-2*oldsol[j]);
                    val = oldsols[j] + (alphaS*dt/(dx*dx))*(oldsols[j+1]-2*oldsols[j] + Lbc)+dt*sources;

                    sols.push_back(val);
                }
                if (j > 0 and j < (NCells-1)){
                    double val;
                    val = oldsols[j] + (alphaS*dt/(dx*dx))*(oldsols[j+1]-2*oldsols[j]+ oldsols[j-1]) + dt*sources;
                    sols.push_back(val);
                }
                if (j==NCells-1){
                    double val;
                    val = oldsols[j] + (alphaf*dt/(dx*dx))*(oldsols[j-1]-2*oldsols[j]+ Rbc) + dt*sources;
                    sols.push_back(val);
                }
                
                //fluid
                if (j==0){
                    double val;
                    val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-Lbc) + (alphaf*dt/(dx*dx))*(oldsolf[j+1]-2*oldsolf[j]+ Lbc)+dt*sourcef;
                    solf.push_back(val);
                }
                if (j > 0 and j < (NCells-1)){
                    double val;
                    val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-oldsolf[j-1]) + (alphaf*dt/(dx*dx))*(oldsolf[j+1]-2*oldsolf[j]+ oldsolf[j-1])+dt*sourcef;
                    solf.push_back(val);
                }
                if (j==NCells-1){
                    double val;
                    val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-oldsolf[j-1]) + (alphaf*dt/(dx*dx))*(oldsolf[j-1]-2*oldsolf[j]+Rbc)+dt*sourcef;
                    solf.push_back(val);
                }
            }
            
            //check if steady state is attained
            
            if (ss_fluid == false){
                double err = 0;
                for(int i = 0; i <NCells;++i){
                    err += 1/dt * abs(solf[i]-oldsolf[i]);
                }
                err = 1.0/NCells*err;
                if(err < errThreshold){
                    std::cerr << "Steady state of fluid solution attained for threshold " << errThreshold << " after " << i << " timesteps." << std::endl;
                    ss_fluid = true;
                }
            }
            
            
            if (ss_solid == false){
                double err = 0;
                for(int i = 0; i < NCells;++i){
                    err += 1/dt * abs(sols[i]-oldsols[i]);
                }
                err = 1.0/NCells*err;
                if(err < errThreshold){
                    std::cerr << "Steady state of solid solution attained for threshold " << errThreshold << " after " << i << " timesteps." << std::endl;
                    ss_solid = true;
                }
            }
            
            //break out of the time step loop when steady state is attained for both solutions
            if(ss_fluid and ss_solid){
                break;
            }
        }
        
        std::vector<double> analyticSol;
        for(int i = 0; i < NCells; ++i){
            analyticSol.push_back(cos(k*dx*i));
        }

        L1F.push_back(L1Error(solf,analyticSol));
        LinfF.push_back(LinfError(solf, analyticSol));
        L1S.push_back(L1Error(sols,analyticSol));
        LinfS.push_back(LinfError(sols, analyticSol));
        
        //update dt
        dt = dt/4;
    }
    exporter.pushOVS(L1F);
    exporter.pushOVS(LinfF);
    exporter.pushOVS(L1S);
    exporter.pushOVS(LinfS);
}

double Simulator::L1Error(std::vector<double> numSol,std::vector<double> analySol){
    if(numSol.size() == analySol.size()){
        double err = 0.0;
        for(int i = 0; i < int(numSol.size()); ++i){
            err += abs((numSol[i]-analySol[i]));
        }
        return (1.0/double(numSol.size()))*err;
    }
    else {
        std::cerr << "L1Error error : the dimensions of the data doesn't match. " << std::endl;
        return 0;
    }
}


double Simulator::LinfError(std::vector<double> numSol,std::vector<double> analySol){
    double maxerror;
    if(numSol.size() == analySol.size()){
        maxerror = abs((numSol[0]-analySol[0]));
        for(int i = 0; i < int(numSol.size()); ++i){
            double val = abs((numSol[i]-analySol[i]));
            if (val > maxerror){
                maxerror = val;
            }
        }
        return maxerror;
    }
    else {
        std::cerr << "LinfError error : the dimensions of the data doesn't match. " << std::endl;
        return 0;
    }
}

//starts the simulation in the charging state
Simulator::Simulator(std::unordered_map<std::string, double> durations,
                     double height,
                     double diameter,
                     int nCells,
                     double T0,
                     int nCycles,
                     int nTimeStepsCycle,
                     Exporter& exporter,
                     double alphaF,
                     double alphaS,
                     double uf):durations(durations),height(height),diameter(diameter),
                    nCells(nCells),T0(T0),nCycles(nCycles),exporter(exporter),nTimeStepsCycle(nTimeStepsCycle),alphaF(alphaF),alphaS(alphaS),uf(uf),state("idledc")
{}
