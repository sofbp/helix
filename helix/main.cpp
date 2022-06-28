#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <sstream>
#include <ostream>
#include <chrono>
#include <map>

#include "peptide.h"
#include "energy_calculator.h"

// 1--> human
// 2--> bacteria

using namespace std;

const double radToDeg = 57.2958;
const double degToRad = 1.0/57.2958;


void crossover(Peptide& pep1, Peptide& pep2){
    //srand (time(NULL));
    int cross_point=rand() %20;

    for (int i = 0; i < cross_point; i++) {
        char temp = pep1.name[i];
        pep1.name[i] = pep2.name[i];
        pep2.name[i]= temp;

    }
}

void mutation(vector<string>& newpopulation,  Peptide pep, Peptide resids){

    //srand (time(NULL));
    int mut_point=rand() %20;

    for (int i=0;i<resids.size();i++) {
        pep.name[mut_point]=resids.name[i];
        string newpeptide=pep.name;
        newpopulation.push_back(newpeptide);
    }

}

void avg_charged(Energy_calculator calc, string input0, string inputCH, vector<vector<double>>& newEh){
        Peptide charged_res;
        charged_res.seq.resize(2);
        calc.Charged.resize(2, vector<vector<double>>(200, vector<double>(2)));
        for (int x=0;x<2;++x) {
            for (int y=0;y<200;++y) {
                for (int z=0;z<2;++z) {
                    calc.Charged[x][y][z]=0;
                }
            }
        }

        calc.load_E_charged(input0, 0);
        calc.load_E_charged(inputCH, 1);

        for (int j=0;j<charged_res.seq.size();j++) {

            charged_res.seq[j].init(calc.Charged[j], 1);

        }
        charged_res.seq[0].E_charged_init(charged_res.seq[1], newEh);

}

void load_charged(Energy_calculator calc, Peptide& resids, map<char, int> res_id, int type){

    if (type==1){

        vector<vector<double>> newEh;
        avg_charged(calc, "K0h", "K+h", newEh);
        resids.seq[res_id['K']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "R0h", "R+h", newEh);
        resids.seq[res_id['R']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "D0h", "D-h", newEh);
        resids.seq[res_id['D']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "E0h", "E-h", newEh);
        resids.seq[res_id['E']].init(newEh, type);
        newEh.clear();

        /*avg_charged(calc, "H0h", "H+h", newEh);
    resids.seq[res_id['H']].init(newEh);
    newEh.clear();*/
    }
    if (type==2){
        vector<vector<double>> newEh;
        avg_charged(calc, "K0b", "K+b", newEh);
        resids.seq[res_id['K']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "R0b", "R+b", newEh);
        resids.seq[res_id['R']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "D0b", "D-b", newEh);
        resids.seq[res_id['D']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "E0b", "E-b", newEh);
        resids.seq[res_id['E']].init(newEh, type);
        newEh.clear();

        /*avg_charged(calc, "H0b", "H+b", newEh);
    resids.seq[res_id['H']].init(newEh);
    newEh.clear();*/
    }
}



int main()
{
    Peptide pep;
    Energy_calculator calc;
    calc.init_grid(20, 200, 2 );

    Peptide resids;
    resids.name="AEKWRCVLIMQFNDSTY";
    resids.seq.resize(resids.name.size());
    map<char, int> res_id{ { 'A', 0 }, { 'E', 1 }, { 'K', 2 },  { 'W', 3 }, { 'R', 4 }, { 'C', 5 }, { 'V', 6 }, { 'L', 7 }, { 'I', 8 },  { 'M', 9 }, { 'Q', 10 }, { 'F', 11 },  { 'N', 12 }, {'D', 13}, {'S', 14}, {'T', 15}, {'Y', 16}};

    for (int j=0;j<resids.name.size();j++) {
        calc.load_EfilesH(resids.name[j], j);
        calc.load_EfilesB(resids.name[j], j);

        resids.seq[j].init(calc.E_H[res_id[resids.name[j]]], 1);
        resids.seq[j].init(calc.E_B[res_id[resids.name[j]]], 2);

    }

    load_charged(calc, resids, res_id, 1);
    load_charged(calc, resids, res_id, 2);

    resids.com_aa("all.gro", res_id);

    pep.load_sequences("input");
    double maxdG=-9999999;
    double max2dG=-9999999;
    Peptide fittest_pep, second_fittest;

    for (int i=0;i<pep.population.size();i++) {
        pep.seq.clear();
        pep.name=pep.population[i];
        pep.seq.resize(20);
        pep.initial_pos("init_pos20", resids, res_id);

        //exit(1);

        //cout<<i<<endl;
        calc.get_energy_total(pep, resids, res_id);

        //calc.print_Emap(pep, resids, res_id);
        //calc.print_Gplot(pep, resids, res_id);
        //exit(0);

        cout<<pep.name<<" "<< pep.energy_h<<" "<<pep.energy_b<<endl;

        // selectivity criteria
        /*if(pep.energy_h>maxdG){
            maxdG=pep.energy_h;
            fittest_pep=pep;
        }
        if(pep.energy_h>max2dG && pep.energy_h < maxdG){
            max2dG=pep.energy_h;
            second_fittest=pep;
        }*/

        //Selectivity criteria with both membranes
        double maximize=pep.energy_h-pep.energy_b;
        if (maximize>maxdG){
            maxdG=maximize;
            fittest_pep=pep;
        }

        if(maximize>max2dG && maximize<maxdG){
            max2dG=maximize;
            second_fittest=pep;
        }
    }


    string peptide =fittest_pep.name;

    //
    // Peptide evolution --> crossover + mutation
    //

    crossover(fittest_pep,second_fittest);

    int rep=0;
    while(rep<100) {

        Peptide newpep;

        mutation(newpep.population,fittest_pep, resids);
        mutation(newpep.population,second_fittest, resids);


        for (int i=0;i<newpep.population.size();i++) {
            newpep.seq.clear();
            newpep.name=newpep.population[i];
            newpep.seq.resize(20);
            newpep.initial_pos("init_pos20", resids, res_id);
            calc.get_energy_total(newpep, resids, res_id);

            //cout<<newpep.name<<" "<< newpep.energy_h<<" "<<newpep.energy_b<<endl;

            // selectivity criteria
            double maximize=newpep.energy_h-newpep.energy_b;
            if (maximize>maxdG){
                maxdG=maximize;
                fittest_pep=newpep;
            }

            if(maximize>max2dG && maximize<maxdG){
                max2dG=maximize;
                second_fittest=newpep;
            }

        }


        string peptide0=fittest_pep.name;

        if(peptide0==peptide){
            rep++;
        }else{
            rep=0;
        }
        cout<<peptide0<<" "<<rep<<" "<<fittest_pep.energy_h<<" "<<fittest_pep.depth_h<<endl;
        peptide=peptide0;
        peptide0.clear();

    }
    calc.print_Emap(fittest_pep, resids, res_id);
    calc.print_Gplot(fittest_pep, resids, res_id);


    return 0;
}
