//
//  Simulator.cpp
//  CFDSim
//
//  Created by Léonard Equer on 02.10.20.
//

#include "Simulator.hpp"
#include <unordered_map>
#include <iostream>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <random>

std::string Simulator::getState(){
    return state;
}

void Simulator::simulate(bool MMS,bool coupled){
    
    //initial temperature
    std::vector<double> sol_fluid(nCells,T0);
    std::vector<double> sol_solid(nCells,T0);
    std::vector<double> oldsol_solid = sol_solid;
    std::vector<double> oldsol_fluid = sol_fluid;
    
    Lbc = 773;
    
    //bool of ss attained, when both true the program stops
    bool ss_fluid = false;
    bool ss_solid = false;
 
    if(MMS){
        Lbc = 1;
        sol_fluid.clear();
        sol_solid.clear();
        for(int nc = 0; nc < nCells-1; ++nc){
            sol_fluid.push_back(cos(this->k_fluid*nc*dx));
            sol_solid.push_back(cos(this->k_solid*nc*dx));
        }
    }
    
    checkStabCond();

    
    for(int i = 1; i < nTimeStepsCycle*nCycles; ++i){
        
        
                
        //export every multiple of pushtimestep solution to a file
        if (i%this->pushTimeStep == 0){
            exporter.pushSolid(sol_solid);
            exporter.pushFluid(sol_fluid);
        }
        
        
        std::vector<double> oldsol_solid = sol_solid;
        std::vector<double> oldsol_fluid = sol_fluid;
        
        sol_solid.clear();
        sol_fluid.clear();
        /*
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
        */
        //get the n+1 solution, pass references to the solution to optimize the speed
        solveDiff(oldsol_solid,oldsol_fluid,sol_solid,sol_fluid,MMS,coupled);
        
        
        //check if steady state is attained
        if(MMS){
            
            //normalize the solid solution
            /*
            double mean = 0;
            for(int i = 0; i < nCells; ++i){
                //calculate average of solid solution to normalize
                mean += (1/double(nCells))*sol_solid[i];
            }
            for(int i = 0; i < nCells; ++i){
                sol_solid[i] -= mean;
            }
            */
            if ((ss_fluid == false) and (i%this->checkSteadyStateTimeStep)){
                double err = 0;
                for(int i = 0; i <nCells;++i){
                    err += 1/dt * abs(sol_fluid[i]-oldsol_fluid[i]);
                }
                err = 1.0/nCells*err;
                if(err < errThreshold){
                    std::cerr << "Steady state of fluid solution attained for threshold " << errThreshold << " after " << i << " timesteps." << std::endl;
                    ss_fluid = true;
                }
            }
            
            if ((ss_solid == false) and (i%this->checkSteadyStateTimeStep)){
                double err = 0;
                for(int i = 0; i < nCells;++i){
                    err += 1/dt * abs(sol_solid[i]-oldsol_solid[i]);
                }
                err = 1.0/nCells*err;
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

    }
    if(MMS){
        std::vector<double> analyticSolFluid;
        std::vector<double> analyticSolSolid;
        

        for(int i = 0; i < nCells; ++i){
            analyticSolFluid.push_back(cos(k_fluid*dx*i));
            analyticSolSolid.push_back(cos(k_solid*dx*i));
        }

        std::cout <<"L1 error fluid: " << L1Error(sol_fluid,analyticSolFluid) << std::endl;
        std::cout << "Linf error fluid: " << LinfError(sol_fluid,analyticSolFluid) << std::endl;
    
        
        std::cout <<"L1 error solid: " << L1Error(sol_solid,analyticSolSolid) << std::endl;
        std::cout << "Linf error solid: " << LinfError(sol_solid,analyticSolSolid) << std::endl;
    }
    
}

void Simulator::solveDiff(const std::vector<double> &oldsol_solid,const  std::vector<double> &oldsol_fluid, std::vector<double> &sol_solid, std::vector<double> &sol_fluid,bool MMS,bool coupled){
    //Method of Manufactured solutions
    double source_solid = 0;
    double source_fluid = 0;
    //puts the right boundary conditions for the MMS
    if(MMS){
        /*
        double xiplu = dx*(1/2);
        double ximin = dx*(-1/2);
        sources = alphaS*k2*(sin(k2*xiplu)-sin(k2*ximin));
        sourcef = alphaF*k*(sin(k*xiplu)-sin(k*ximin)) + uf*(cos(k*xiplu) - cos(k*ximin));
        */
        double xi = 0;
        source_solid = alphaS*k_solid*k_solid*cos(k_solid*xi);
        source_fluid = alphaF*k_fluid*k_fluid*cos(k_fluid*xi) - uf*k_fluid*sin(k_fluid*xi);
        if(coupled){
            source_solid += hvs*(cos(k_solid*xi)-cos(k_fluid*xi));
            source_fluid += hvf*(cos(k_fluid*xi)-cos(k_solid*xi));
        }
    }

    //left boundary values
    sol_solid.push_back(oldsol_solid[0] + (alphaS*dt/(dx*dx))*(oldsol_solid[1]-oldsol_solid[0])+(dt)*source_solid);
    sol_fluid.push_back(oldsol_fluid[0] - (uf*dt/dx)*(oldsol_fluid[0]-Lbc) + (alphaF*dt/(dx*dx))*(oldsol_fluid[1]-oldsol_fluid[0])+(dt)*source_fluid);
    

    //over all space for each timestep
    for(int j = 1; j<=nCells-2;++j){
        
        //calculate the source term in the MMS method
        if(MMS){
            /*
            double xiplu = dx*(j+1/2);
            double ximin = dx*(j-1/2);
            sources = alphaS*k2*(sin(k2*xiplu)-sin(k2*ximin));
            sourcef = alphaF*k*(sin(k*xiplu)-sin(k*ximin)) + uf*(cos(k*xiplu) - cos(k*ximin));
            */
            double xi = dx*(j);
            source_solid = alphaS*k_solid*k_solid*cos(k_solid*xi);
            source_fluid = alphaF*k_fluid*k_fluid*cos(k_fluid*xi) - uf*k_fluid*sin(k_fluid*xi);
            if(coupled){
                source_solid += hvs*(cos(k_solid*xi)-cos(k_fluid*xi));
                source_fluid += hvf*(cos(k_fluid*xi)-cos(k_solid*xi));
            }
        }
        
        //solid
        sol_solid.push_back(oldsol_solid[j] + (alphaS*dt/(dx*dx))*(oldsol_solid[j+1]-2*oldsol_solid[j]+ oldsol_solid[j-1]) + (dt)*source_solid);
 
        //fluid
        sol_fluid.push_back(oldsol_fluid[j] - (uf*dt/dx)*(oldsol_fluid[j]-oldsol_fluid[j-1]) + (alphaF*dt/(dx*dx))*(oldsol_fluid[j+1]-2*oldsol_fluid[j]+ oldsol_fluid[j-1])+(dt)*source_fluid);
        
        if(coupled){
            //solve the coupled equation one step behind to optimize memory
            Eigen::Matrix2d A;
            A << (1+hvf*dt), (-hvf*dt), (-hvs*dt), (1+hvs*dt);
            Eigen::Vector2d b(sol_fluid[j-1],sol_solid[j-1]);
            Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
            
            //puts the coupled solution back in
            sol_solid[j-1] = x(1);
            sol_fluid[j-1] = x(0);
        }
    }
    
    //MMS last value
    if(MMS){
        /*
        double xiplu = dx*(nCells - 1 + 1/2);
        double ximin = dx*(nCells - 1 - 1/2);
        sources = alphaS*k2*(sin(k2*xiplu)-sin(k2*ximin));
        sourcef = alphaF*k*k*(sin(k*xiplu)-sin(k*ximin)) + uf*k*(cos(k*xiplu)- cos(k*ximin));
         */
        
        double xi = dx*(nCells-1);
        source_solid = alphaS*k_solid*k_solid*cos(k_solid*xi);
        source_fluid = alphaF*k_fluid*k_fluid*cos(k_fluid*xi) - uf*k_fluid*sin(k_fluid*xi);
        if(coupled){
            source_solid += hvs*(cos(k_solid*xi)-cos(k_fluid*xi));
            source_fluid += hvf*(cos(k_fluid*xi)-cos(k_solid*xi));
        }
    }
    
    //values on the right boundary
    sol_solid.push_back(oldsol_solid[nCells-1] + (alphaS*dt/(dx*dx))*(oldsol_solid[nCells-2]-oldsol_solid[nCells-1]) + (dt)*source_solid);
    
    sol_fluid.push_back(oldsol_fluid[nCells-1] - (uf*dt/dx)*(oldsol_fluid[nCells-1]-oldsol_fluid[nCells-2]) + (alphaF*dt/(dx*dx))*(oldsol_fluid[nCells-2]-oldsol_fluid[nCells-1])+(dt)*source_fluid);
    
    if(coupled){
        //solve the coupled equation one step behind to optimize memory
        Eigen::Matrix2d A;
        A << (1+hvf*dt), (-hvf*dt), (-hvs*dt), (1+hvs*dt);
        Eigen::Vector2d b(sol_fluid[nCells-2],sol_solid[nCells-2]);
        Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
;
        
        //puts the coupled solution back in
        sol_solid[nCells-2] = x(1);
        sol_fluid[nCells-2] = x(0);
        
        b(0) = sol_fluid[nCells-1];
        b(1) = sol_solid[nCells-1];
        x = A.colPivHouseholderQr().solve(b);
        sol_solid[nCells-1] = x(1);
        sol_fluid[nCells-1] = x(0);
    }
    
    
}

void Simulator::checkStabCond(){
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
    }
    
    if(sigma + 2*df > 1){
        std::cerr << "The method might not be stable (upper bound)"  << std::endl;
    }
}



void Simulator::OVSNonCoupledDiff(double Pe, int n){
    this->dt = 20;
    this -> n_fluid = n;
    this -> n_solid = n;
    
    
    //std::default_random_engine generator;
    //std::uniform_real_distribution<double> dist(-0.1,0.1);
    
    std::vector<double> L1F;
    std::vector<double> LinfF;
    std::vector<double> L1S;
    std::vector<double> LinfS;
    
    int NCellsi = 8;
    
    //puts the right boundary conditions for the MMS
    this->Lbc = 1;
    
    
    for(int i = 0;i < 8; ++i){
        
        this->nCells = NCellsi*pow(2,i);
        
        std::cout << "Number of cells : " << nCells << std::endl;
        //param
        this->dx = height/double(nCells);


        std::cout << "dt : " << dt << std::endl;
        std::cout << "dx : " << dx << std::endl;
        
        //this->alphaF = 0.03 * dx *dx /dt;
        this->uf = (alphaF*Pe)/(dx*nCells);
        //this->alphaS = 0.03 * dx *dx /dt;

 
        //bool of ss attained, when both true the program stops
        bool ss_fluid = false;
        bool ss_solid = false;
        
        //to have these variables as global
        std::vector<double> solf;
        std::vector<double> sols;
        
        //initial temperature
        for(int nc = 0; nc < nCells-1; ++nc){
            solf.push_back(cos(this->k_fluid*nc*dx));
            sols.push_back(cos(this->k_solid*nc*dx));
            //solf.push_back(cos(this->k_fluid*nc*dx)+dist(generator));
            //sols.push_back(cos(this->k_solid*nc*dx)+dist(generator));
        }
        
        std::vector<double> oldsols = sols;
        std::vector<double> oldsolf = solf;
    
        //check stability conditions
        checkStabCond();
        
        //timestep simulation of the program
        for(int i = 1; i <= nTimeStepsCycle; ++i){
            
            std::vector<double> oldsols = sols;
            std::vector<double> oldsolf = solf;

            sols.clear();
            solf.clear();
            
            //get the n+1 solution, pass references to the solution to optimize the speed
            solveDiff(oldsols,oldsolf,sols,solf,true,false);
            /*
            double mean = 0;
            for(int i = 0; i < nCells; ++i){
                //calculate average of solid solution to normalize
                mean += (1/double(nCells))*sols[i];
            }
            for(int i = 0; i < nCells; ++i){
                sols[i] -= mean;
            }
            */
            //check if steady state is attained
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
            
            //break out of the time step loop when steady state is attained for both solutions
            if(ss_fluid and ss_solid){
                break;
            }
        }
        
        std::vector<double> analyticSolFluid;
        std::vector<double> analyticSolSolid;

        for(int i = 0; i < nCells; ++i){
            analyticSolFluid.push_back(cos(k_fluid*dx*i));
            analyticSolSolid.push_back(cos(k_solid*dx*i));
        }

        L1F.push_back(L1Error(solf,analyticSolFluid));
        LinfF.push_back(LinfError(solf, analyticSolFluid));
        
        L1S.push_back(L1Error(sols,analyticSolSolid));
        LinfS.push_back(LinfError(sols, analyticSolSolid));
        
        //update dt
        dt = dt/2;
        
        
        std::cout << "L1F : " << L1F[i] << std::endl;
        std::cout << "LinfF : " <<LinfF[i] << std::endl;
        std::cout << "L1S : " <<L1S[i] << std::endl;
        std::cout << "LinfS : " <<LinfS[i] << std::endl;
    }


    exporter.pushOVS(L1F);
    exporter.pushOVS(LinfF);
    
    exporter.pushOVS(L1S);
    exporter.pushOVS(LinfS);
    
}

double Simulator::L1Error(const std::vector<double> &numSol,const std::vector<double> &analySol){
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


double Simulator::LinfError(const std::vector<double> &numSol,const std::vector<double> &analySol){
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
                     double uf,
                     double hvs,
                     double hvf):durations(durations),height(height),diameter(diameter),
                    nCells(nCells),T0(T0),nCycles(nCycles),exporter(exporter),nTimeStepsCycle(nTimeStepsCycle),alphaF(alphaF),alphaS(alphaS),uf(uf),state("idledc"),hvs(hvs),hvf(hvf)
{
    
    //seting up some constants
    dx = height/double(nCells);
    std::cout << "dt : " << dt << std::endl;
    std::cout << "dx : " << dx << std::endl;
    k_fluid = 2*M_PI*n_fluid/this->height;
    k_solid = 2*M_PI*n_solid/this->height;
}

/*
 ******************************************************************************
 Deprecated, for reference
 *******************************************************************************
*/

void Simulator::solveNonCoupledDiff(bool MMS){

    //boundary conditions
    double Lbc = 500;

    //bool of ss attained, when both true the program stops
    bool ss_fluid = false;
    bool ss_solid = false;
    
    //initial temperature
    std::vector<double> solf(nCells,T0);
    std::vector<double> sols(nCells,T0);
    if(MMS){
        solf.clear();
        sols.clear();
        for(int nc = 0; nc < nCells - 1; ++nc){
            solf.push_back(cos(this->k_fluid*nc*dx));
            sols.push_back(cos(this->k_fluid*nc*dx));
        }
    }
    
    //to have these variables as global
    std::vector<double> oldsols = sols;
    std::vector<double> oldsolf = solf;
    
    
    //Method of Manufactured solutions
    
    double sources = 0;
    double sourcef = 0;
    //puts the right boundary conditions for the MMS
    if(MMS){
        Lbc = 1;
    }
    
    //check stability conditions
    checkStabCond();
    
    
    
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

        for(int j = 0; j<=nCells-1;++j){
            
            //calculate the source term in the MMS method
            if(MMS){
                double xi = dx*j;
                sources = alphaS*k_fluid*k_fluid*cos(k_fluid*xi);
                sourcef = alphaF*k_fluid*k_fluid*cos(k_fluid*xi) - uf*k_fluid*sin(k_fluid*xi);
            }
            
            //solid
            if (j==0){
                double val;
                //val = oldsol[j] + (alphaS*dt/(dx*dx))*(oldsol[j+1]-2*oldsol[j]);
                val = oldsols[j] + (alphaS*dt/(dx*dx))*(oldsols[j+1]-oldsols[j])+dt*sources;

                sols.push_back(val);
            }
            if (j > 0 and j < (nCells-1)){
                double val;
                val = oldsols[j] + (alphaS*dt/(dx*dx))*(oldsols[j+1]-2*oldsols[j]+ oldsols[j-1]) + dt*sources;
                sols.push_back(val);
            }
            if (j==nCells-1){
                double val;
                val = oldsols[j] + (alphaF*dt/(dx*dx))*(oldsols[j-1]-oldsols[j]) + dt*sources;
                sols.push_back(val);
            }
            
            //fluid
            if (j==0){
                double val;
                val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-Lbc) + (alphaF*dt/(dx*dx))*(oldsolf[j+1]-oldsolf[j])+dt*sourcef;
                solf.push_back(val);
            }
            if (j > 0 and j < (nCells-1)){
                double val;
                val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-oldsolf[j-1]) + (alphaF*dt/(dx*dx))*(oldsolf[j+1]-2*oldsolf[j]+ oldsolf[j-1])+dt*sourcef;
                solf.push_back(val);
            }
            if (j==nCells-1){
                double val;
                val = oldsolf[j] - (uf*dt/dx)*(oldsolf[j]-oldsolf[j-1]) + (alphaF*dt/(dx*dx))*(oldsolf[j-1]-oldsolf[j])+dt*sourcef;
                solf.push_back(val);
            }
        }
        
        //check if steady state is attained
        if(MMS){
            
            if ((ss_fluid == false) and (i%this->checkSteadyStateTimeStep)){
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
            
            if ((ss_solid == false) and (i%this->checkSteadyStateTimeStep)){
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
            //break out of the time step loop when steady state is attained for both solutions
            if(ss_fluid and ss_solid){
                break;
            }
        }
        
    }
    if(MMS){
        
        std::vector<double> analyticSol;
        for(int i = 0; i < nCells; ++i){
            analyticSol.push_back(cos(k_fluid*dx*i));
        }

        std::cout <<"L1 error : " << L1Error(solf,analyticSol) << std::endl;
        std::cout << "Linf error : " << LinfError(solf,analyticSol) << std::endl;
    }
}
