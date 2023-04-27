#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#include "fileutils.h"
#include "stringutils.h"
#include "gnuplot.h"

using namespace std;


// ==============================================================================================================================
// Define a function to calculate the energy of the grid using the Ising model Hamiltonian
double getEnergy(const vector<vector<int>>& grid) {
    // printf("getEnergy");
    int N = grid.size();
    double energy = 0.0;
    double J = 0.5;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int spin = grid[i][j];
            int neighbors = grid[(i-1+N)%N][j] + grid[(i+1)%N][j] + grid[i][(j-1+N)%N] + grid[i][(j+1)%N];
            energy += -J*spin * neighbors;
        }
    }
    return energy / 2.0;
}

// ==============================================================================================================================
// Define a function to calculate the total magnetization of the grid
int getMagnetization(const vector<vector<int>>& grid) {
    int N = grid.size();
    int magnetization = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            magnetization += grid[i][j];
        }
    }
    return magnetization;
}

    

// ==============================================================================================================================
// Define the main function that implements the Metropolis algorithm
vector<vector<int>> ising_model(int N, double temperature, int num_steps,fHandle fE,string sE) {
    
    

    // Step 1: Initialize a random grid
    vector<vector<int>> grid(N, vector<int>(N));
    
    
    random_device rd;
    mt19937 gen(rd()); // random number generator
    uniform_int_distribution<> dis(0, 1); //creates a uniform distribution between 0 and 1
    
    
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            grid[i][j] = 2*dis(gen) - 1; // fill in random values 
        }
    }


    // Step 2: Define variables to keep track of energy and magnetization
    double energy = getEnergy(grid);
    int magnetization = getMagnetization(grid);
    double e,m; 

    // Step 3: Define the Metropolis algorithm
    for (int step = 0; step < num_steps; step++) {

        // Step 4: Choose a random spin
        uniform_int_distribution<> dis2(0, N-1);


        
        int i = dis2(gen);
        int j = dis2(gen);
        int localSpin = grid[i][j];

        // int iup = i+1;
        // int idown = i-1;
        // int jup = j+1;
        // int jdown = j-1;
        
        // Calculate the change in energy if the spin is flipped
        int neighbors = grid[(i-1+N)%N][j] + grid[(i+1)%N][j] + grid[i][(j-1+N)%N] + grid[i][(j+1)%N];
        //int neighbors = grid[(idown)][j] + grid[iup][j] + grid[i][jdown] + grid[i][jup];
        double delta_energy = 2.0 * localSpin * neighbors; // comes from the change in energy when flipping from -1 to 1.

        // Accept or reject the flip based on the Metropolis criterion
		//if change in energy is negative, or if the energy is greater than some random number, change it

        
        if (delta_energy < 0.0) {
            if(dis(gen) < exp(-delta_energy / temperature)){
                localSpin *= -1; //flip spin
                grid[i][j] = localSpin; //update it
                energy += delta_energy; //add to the total energy the change in energy
                magnetization += 2 * localSpin; //add to the total magnetization the changed spin
                //printf("energy : %f \t magnetization: %d \n", energy, magnetization);
            }
        } 
        // sE = FloatToStr(step) + "\t"+FloatToStr(energy) + "\t"+FloatToStr(magnetization) + "\n";
        // FileWrite(fE,sE.c_str(),sE.length());
    }
    return grid;
}

// ==============================================================================================================================
// Example usage of the ising_model function
int main() {
    int N = 10;
    int N2 = N*N;
    
    double temperature = 3.0;
    int num_steps = 1000000;
    gnuplot *gp = new gnuplot();
    
    string sC;
    fHandle fC;
    fC = FileCreate("config.txt");
    
    string sE;
    fHandle fE;
    fE = FileCreate("energy.txt");

    int spin;
    
    vector<vector<int>> grid = ising_model(N, temperature, num_steps,fE,sE); //this reaches steady state? 


    //writing the equilibrated spin configuration to file for graphing 
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            sC = IntToStr(grid[i][j]) + "\t";
            FileWrite(fC,sC.c_str(),sC.length());
        }
        sC = "\n";
        FileWrite(fC,sC.c_str(),sC.length());
    }
    FileClose(fC);

    double e,e2;
    int m;
    double avgE; //<E>
    double avgM; //<M>
    double cH; //specific heat capacity
    double X;//suseptibility

    double overN2 = 1/((double) N2);
    double overTemp;
    string s;
    fHandle f;
    f = FileCreate("ising.txt");
    for(int t = 1; t < 101; t++) { //averages vs temperature plots
        temperature = 0.03*t;
        overTemp = 1/temperature;
        printf("Current temperature is: %f\n",temperature);

        //here is where I would want to iterate through and calclate / write to file all the measurements
        grid  = ising_model(N, temperature, num_steps,fE, sE); 
        e = getEnergy(grid); // energy for the config
        e2 = e*e;
        m = getMagnetization(grid); //magnetization for the config
        

        avgE = e * overN2;   
        // printf("e: %f, overN2: %f, avgE: %f \n",e, overN2,avgE);      
        avgM = ((double) m )* overN2;
        cH = (e*e*overN2 - (avgE * avgE)) * overTemp;
        X = (m*m*overN2 - (avgM * avgM))*overTemp;
        s = FloatToStr(temperature) + "\t"+FloatToStr(avgE) + "\t"+FloatToStr(avgM) + "\t"+FloatToStr(cH) + "\t"+FloatToStr(X) + "\n";
        FileWrite(f,s.c_str(),s.length());
    }
    FileClose(fE);
    FileClose(f);    

    gp->plotfile("ising.txt","u 1:2 w l t '<E>'");
    gp->plotfile("ising.txt","u 1:3 w l t '<M>'");
    gp->plotfile("ising.txt","u 1:4 w l t 'c_H'");
    gp->plotfile("ising.txt","u 1:5 w l t 'X'");
    gp->show();

    // gp->plotfile("ising.txt","u 1:2 w l t '<E>'");
    // gp->replotfile("ising.txt","u 1:3 w l t '<M>'");
    // gp->show();

    delete gp;

    return 0;
}