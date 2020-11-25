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
#include <string>
#include <fstream>

std::string Simulator::getState(){
    return state;
}

void Simulator::simulate(bool MMS,bool coupled){

    //initial temperature
    std::vector<double> sol_fluid(nCells,T0);
    std::vector<double> sol_solid(nCells,T0);
    std::vector<double> oldsol_solid = sol_solid;
    std::vector<double> oldsol_fluid = sol_fluid;


    sim_uf = this->uf;
    
    //bool of ss attained, when both true the program stops
    bool ss_fluid = false;
    bool ss_solid = false;

    if(MMS){
        Lbc = 1;
        sol_fluid.clear();
        sol_solid.clear();
        for(int nc = 0; nc < nCells; ++nc){
            sol_fluid.push_back(cos(this->k_fluid*nc*dx));
            sol_solid.push_back(cos(this->k_solid*nc*dx));
        }
    }

    checkStabCond();


    for(int i = 1; i < nTimeStepsCycle*nCycles; ++i){
        //export every multiple of pushtimestep solution to a file
        if ((i%this->pushTimeStep == 0) or (i == nTimeStepsCycle*nCycles-1)){
            exporter.pushSolid(sol_solid);
            exporter.pushFluid(sol_fluid);
        }

        std::vector<double> oldsol_solid = sol_solid;
        std::vector<double> oldsol_fluid = sol_fluid;

        sol_solid.clear();
        sol_fluid.clear();

        int total_time_cycle = durations["discharging"]+durations["idlecd"]+durations["charging"] + durations["idledc"];
        double val = i%nTimeStepsCycle;
        //updating the state of the simulation
        if (val == nTimeStepsCycle*durations["charging"]/total_time_cycle){
            this->state = "idlecd";
            sim_uf = 0;
            if(!MMS){
                std::cout << "Temperature difference between beginning and end of charging : " << oldsol_fluid.back()-T0 << std::endl;
            }
        }
        if (val == nTimeStepsCycle*(durations["idlecd"]+durations["charging"])/total_time_cycle){
            this->state = "discharging";
            sim_uf = -this->uf ;
            ThermalEnergyBefDischar = computeStoredThermalEnergy(sol_solid, sol_fluid);

        }
        if (val == nTimeStepsCycle*(durations["discharging"]+durations["idlecd"]+durations["charging"])/total_time_cycle){
            this->state = "idledc";
            sim_uf = 0;
            
            if(!MMS){
                //capacity factors
                double capa = (ThermalEnergyBefDischar-ThermalEnergyBefChar)/maximumThermalEnergy;
                capacityFactors.push_back(capa);
                std::cout << "Capacity factor for this cycle : ";
                std::cout << capa << std::endl;
                
                //exergy efficiencies
                double exEff = (fluxdout-fluxdin)/(fluxcin-fluxcout);
                exergyEfficiencies.push_back(exEff);
                std::cout << "Exergy efficiency for this cycle : " << exEff << std::endl;
                fluxdout = 0;
                fluxcin = 0;
                fluxdin = 0;
                fluxcout = 0;
                
            }
            
        }
        if (val == 0){
            this->state = "charging";
            sim_uf = this->uf;
            ThermalEnergyBefChar = computeStoredThermalEnergy(sol_solid, sol_fluid);
        
        }

        exporter.exportState(state);

        //makes sure we are in charging state when doing MMS
        if(MMS){
            this->state = "charging";
            sim_uf = this->uf;
        }

        //get the n+1 solution, pass references to the solution to optimize the speed
        solveDiff(oldsol_solid,oldsol_fluid,sol_solid,sol_fluid,MMS,coupled);
        
        //check if steady state is attained
        if(MMS){

            if ((ss_fluid == false) and (i%this->checkSteadyStateTimeStep)){
                double err = 0;
                err = L1Error(sol_fluid, oldsol_fluid);
                err = err/dt;
                if(err < errThreshold){
                    std::cerr << "Steady state of fluid solution attained for threshold " << errThreshold << " after " << i << " timesteps." << std::endl;
                    ss_fluid = true;
                }
            }

            if ((ss_solid == false) and (i%this->checkSteadyStateTimeStep)){
                double err = 0;
                err = L1Error(sol_solid, oldsol_solid);
                err = err/dt;
                if(err < errThreshold){
                    std::cerr << "Steady state of solid solution attained for threshold " << errThreshold << " after " << i << " timesteps." << std::endl;
                    ss_solid = true;
                }
            }
            //break out of the time step loop when steady state is attained for both solutions
            if(ss_fluid and ss_solid){
                break;
            }
        }else{
            //get the next value of the integral for the exergy flux calculation
            computeNextExergyFlux(oldsol_fluid,sol_fluid);
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
        source_fluid = alphaF*k_fluid*k_fluid*cos(k_fluid*xi) - sim_uf*k_fluid*sin(k_fluid*xi);
        if(coupled){
            source_solid += hvs*(cos(k_solid*xi)-cos(k_fluid*xi));
            source_fluid += hvf*(cos(k_fluid*xi)-cos(k_solid*xi));
        }
    }
    if(state == "charging" or state == "idlecd"){
        //left boundary values
        sol_solid.push_back(oldsol_solid[0] + (alphaS*dt/(dx*dx))*(oldsol_solid[1]-oldsol_solid[0])+(dt)*source_solid);
        sol_fluid.push_back(oldsol_fluid[0] - (sim_uf*dt/dx)*(oldsol_fluid[0]-Lbc) + (alphaF*dt/(dx*dx))*(oldsol_fluid[1]-oldsol_fluid[0])+(dt)*source_fluid);
    }
    if(state == "discharging" or state == "idledc"){
        //left boundary values
        sol_solid.push_back(oldsol_solid[0] + (alphaS*dt/(dx*dx))*(oldsol_solid[1]-oldsol_solid[0])+(dt)*source_solid);
        sol_fluid.push_back(oldsol_fluid[0] - (sim_uf*dt/dx)*(oldsol_fluid[1]-oldsol_fluid[0]) + (alphaF*dt/(dx*dx))*(oldsol_fluid[1]-oldsol_fluid[0])+(dt)*source_fluid);
    }


    //over all space for each timestep
    for(int j = 1; j<=nCells-2;++j){

        //calculate the source term in the MMS method
        if(MMS){
            /*
            double xiplu = dx*(j+1/2);
            double ximin = dx*(j-1/2);
            sources = alphaS*k2*(sin(k2*xiplu)-sin(k2*ximin));
            sourcef = alphaF*k*(sin(k*xiplu)-sin(k*ximin)) + sim_uf*(cos(k*xiplu) - cos(k*ximin));
            */
            double xi = dx*(j);
            source_solid = alphaS*k_solid*k_solid*cos(k_solid*xi);
            source_fluid = alphaF*k_fluid*k_fluid*cos(k_fluid*xi) - sim_uf*k_fluid*sin(k_fluid*xi);
            if(coupled){
                source_solid += hvs*(cos(k_solid*xi)-cos(k_fluid*xi));
                source_fluid += hvf*(cos(k_fluid*xi)-cos(k_solid*xi));
            }
        }

        if(state == "charging" or state == "idlecd"){
            //solid
            sol_solid.push_back(oldsol_solid[j] + (alphaS*dt/(dx*dx))*(oldsol_solid[j+1] - 2*oldsol_solid[j] + oldsol_solid[j-1]) + (dt)*source_solid);                

            //fluid
            sol_fluid.push_back(oldsol_fluid[j] - (sim_uf*dt/dx)*(oldsol_fluid[j]-oldsol_fluid[j-1]) + (alphaF*dt/(dx*dx))*(oldsol_fluid[j+1]-2*oldsol_fluid[j] + oldsol_fluid[j-1])+(dt)*source_fluid);
        }
        if(state == "discharging" or state == "idledc"){
            //solid
            sol_solid.push_back(oldsol_solid[j] + (alphaS*dt/(dx*dx))*(oldsol_solid[j+1]-2*oldsol_solid[j]+ oldsol_solid[j-1]) + (dt)*source_solid);

            //fluid
            sol_fluid.push_back(oldsol_fluid[j] - (sim_uf*dt/dx)*(oldsol_fluid[j+1]-oldsol_fluid[j]) + (alphaF*dt/(dx*dx))*(oldsol_fluid[j+1]-2*oldsol_fluid[j]+ oldsol_fluid[j-1])+(dt)*source_fluid);
        }


        if(coupled){
            //solve the coupled equation one step behind to optimize memory
            
            double new_Tf = (sol_fluid[j-1] + hvf*dt*sol_solid[j-1])/(1+hvf*dt);
            double new_Ts = (sol_solid[j-1] + hvs*dt*sol_fluid[j-1]/(1+hvf*dt))/(1+hvs*dt-(hvs*hvf*dt*dt/(1+hvf*dt)));


            //puts the coupled solution back in
            sol_solid[j-1] = new_Ts;
            sol_fluid[j-1] = new_Tf;
        }
    }

    //MMS last value
    if(MMS){
        /*
        double xiplu = dx*(nCells - 1 + 1/2);
        double ximin = dx*(nCells - 1 - 1/2);
        sources = alphaS*k2*(sin(k2*xiplu)-sin(k2*ximin));
        sourcef = alphaF*k*k*(sin(k*xiplu)-sin(k*ximin)) + sim_uf*k*(cos(k*xiplu)- cos(k*ximin));
         */

        double xi = dx*(nCells-1);
        source_solid = alphaS*k_solid*k_solid*cos(k_solid*xi);
        source_fluid = alphaF*k_fluid*k_fluid*cos(k_fluid*xi) - sim_uf*k_fluid*sin(k_fluid*xi);
        if(coupled){
            source_solid += hvs*(cos(k_solid*xi)-cos(k_fluid*xi));
            source_fluid += hvf*(cos(k_fluid*xi)-cos(k_solid*xi));
        }
    }

    if(state == "charging" or state == "idlecd"){
        //values on the right boundary
        sol_solid.push_back(oldsol_solid[nCells-1] + (alphaS*dt/(dx*dx))*(oldsol_solid[nCells-2]-oldsol_solid[nCells-1]) + (dt)*source_solid);

        sol_fluid.push_back(oldsol_fluid[nCells-1] - (sim_uf*dt/dx)*(oldsol_fluid[nCells-1]-oldsol_fluid[nCells-2]) + (alphaF*dt/(dx*dx))*(oldsol_fluid[nCells-2]-oldsol_fluid[nCells-1])+(dt)*source_fluid);
    }
    if(state == "discharging" or state == "idledc"){
        //values on the right boundary
        sol_solid.push_back(oldsol_solid[nCells-1] + (alphaS*dt/(dx*dx))*(oldsol_solid[nCells-2]-oldsol_solid[nCells-1]) + (dt)*source_solid);

        sol_fluid.push_back(oldsol_fluid[nCells-1] - (sim_uf*dt/dx)*(Rbc-oldsol_fluid[nCells-1]) + (alphaF*dt/(dx*dx))*(oldsol_fluid[nCells-2]-oldsol_fluid[nCells-1])+(dt)*source_fluid);
    }


    if(coupled){
        //solve the coupled equation one step behind to optimize memory
        
        double new_Tf = (sol_fluid[nCells-2] + hvf*dt*sol_solid[nCells-2])/(1+hvf*dt);
        double new_Ts = (sol_solid[nCells-2] + hvs*dt*sol_fluid[nCells-2]/(1+hvf*dt))/(1+hvs*dt-(hvs*hvf*dt*dt/(1+hvf*dt)));

        //puts the coupled solution back in
        sol_solid[nCells-2] = new_Ts;
        sol_fluid[nCells-2] = new_Tf;
        
        new_Tf = (sol_fluid[nCells-1] + hvf*dt*sol_solid[nCells-1])/(1+hvf*dt);
        new_Ts = (sol_solid[nCells-1] + hvs*dt*sol_fluid[nCells-1]/(1+hvf*dt))/(1+hvs*dt-(hvs*hvf*dt*dt/(1+hvf*dt)));

        sol_solid[nCells-1] = new_Ts;
        sol_fluid[nCells-1] = new_Tf;
    }


}

void Simulator::checkStabCond(){
    
    std::cout << std::endl;
    std::cout  << "CFL and diffusion numbers : " << std::endl;
    //check stability conditions
    double sigma = sim_uf*dt/dx;
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
    std::cout << std::endl;
}



void Simulator::OVS(double Pe, int n, bool coupled){
    this->dt = 5;
    this -> n_fluid = n;
    this -> n_solid = n;
    k_fluid = 2*M_PI*n_fluid/this->height;
    k_solid = 2*M_PI*n_solid/this->height;
    if(coupled){
        this -> n_fluid = 2;
        this -> n_solid = 1;
        
        std::cout << "Coupled equations OVS with Pe = " << Pe << std::endl;
    }else{
        std::cout << "Non coupled equations OVS with Pe = " << Pe << " and n = " << n <<  std::endl;
    }

    std::vector<double> L1F;
    std::vector<double> LinfF;
    std::vector<double> L1S;
    std::vector<double> LinfS;

    int NCellsi = 16;

    //puts the right boundary conditions for the MMS
    this->Lbc = 1;
    this->sim_uf = uf;

    for(int i = 0;i < 5; ++i){

        this->nCells = NCellsi*pow(2,i);

        std::cout << "Number of cells : " << nCells << std::endl;
        //param
        this->dx = height/double(nCells);


        std::cout << "dt : " << dt << std::endl;
        std::cout << "dx : " << dx << std::endl;

        //this->alphaF = 0.03 * dx *dx /dt;
        this->sim_uf = (alphaF*Pe)/(this->height);
        //this->alphaS = 0.03 * dx *dx /dt;


        //bool of ss attained, when both true the program stops
        bool ss_fluid = false;
        bool ss_solid = false;

        //to have these variables as global
        std::vector<double> solf;
        std::vector<double> sols;

        //initial temperature
        for(int nc = 0; nc < nCells; ++nc){
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
            solveDiff(oldsols,oldsolf,sols,solf,true,coupled);
            
            //check if steady state is attained
            if (ss_fluid == false){
                double err = 0;
                err = L1Error(solf, oldsolf);
                err = err/dt;
                if(err < errThreshold){
                    std::cerr << "Steady state of fluid solution attained for threshold " << errThreshold << " after " << i << " timesteps." << std::endl;
                    ss_fluid = true;
                }
            }


            if (ss_solid == false){
                double err = 0;
                err = L1Error(sols, oldsols);
                err = err/dt;
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
        std::cerr << numSol.size() <<" " << analySol.size() << std::endl;
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

double Simulator::computeStoredThermalEnergy(std::vector<double> &sols, std::vector<double> &solf){
    //trapezoidal formula
    
    double integral_fluid = 0;
    double integral_solid = 0;

    for(int i = 0; i < nCells-1;++i){
        integral_fluid += (((i+1)*dx-i*dx)/2)*(solf[i+1]+solf[i]);
        integral_solid += (((i+1)*dx-i*dx)/2)*(sols[i+1]+sols[i]);
    }
    
    double Qt =(M_PI/4)*pow(diameter,2)*(0.4*1835.6*1511.8*integral_fluid+(1-0.4)*2600*900*integral_solid);
    return Qt;
}

void Simulator::computeNextExergyFlux(std::vector<double> &oldsolf,std::vector<double> &solf){
    //uses trapezoidal at every timestep during the charging and discharging phases to update the integral values
    double prefactor = 10.0*1511.8;
    double T = solf[0];
    double oldT = oldsolf[0];
    if(state == "discharging"){
        fluxdin += prefactor*dt*0.5*((T-T0-T0*std::log(T/T0))+(oldT-T0-T0*std::log(oldT/T0)));
    }else if (state == "charging"){
        fluxcin += prefactor*dt*0.5*((T-T0-T0*std::log(T/T0))+(oldT-T0-T0*std::log(oldT/T0)));
    }
    T = solf.back();
    oldT = oldsolf.back();

    if(state == "discharging"){
        fluxdout += prefactor*dt*0.5*((T-T0-T0*std::log(T/T0))+(oldT-T0-T0*std::log(oldT/T0)));
    }else if (state == "charging"){
        fluxcout += prefactor*dt*0.5*((T-T0-T0*std::log(T/T0))+(oldT-T0-T0*std::log(oldT/T0)));
    }
    
}


//starts the simulation in the charging state
Simulator::Simulator(std::unordered_map<std::string, int> durations,
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
                    nCells(nCells),T0(T0),nCycles(nCycles),exporter(exporter),nTimeStepsCycle(nTimeStepsCycle),alphaF(alphaF),alphaS(alphaS),uf(uf),state("charging"),hvs(hvs),hvf(hvf)
{

    //seting up some constants
    dx = height/double(nCells);
    dt = (double(durations["idledc"])+double(durations["discharging"])+double(durations["idlecd"])+double(durations["charging"]))/double(nTimeStepsCycle);

    std::cout << "dt : " << dt << std::endl;
    std::cout << "dx : " << dx << std::endl;
    k_fluid = 2*M_PI*n_fluid/this->height;
    k_solid = 2*M_PI*n_solid/this->height;
    
    Lbc = 873;
    Rbc = 293;
  
    //setting up the maximum thermal energy
    double prefactor = 0.4*1835.6*1511.8+(1-0.4)*2600*900;
    this->maximumThermalEnergy = prefactor*(M_PI/4)*pow(this->diameter,2)*this->height*(this->Lbc-this->Rbc);
    


}
