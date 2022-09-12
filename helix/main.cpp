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


void mutation(vector<string>& newpopulation,  string name, Peptide resids){

    srand (time(NULL));
    int mut_point=rand() %20;
    int before=mut_point-4;
    int before1=mut_point-3;
    int after=mut_point+4;
    int after1=mut_point+3;

    if(((name[before] == 'K' || name[before1] == 'K') && (name[after] == 'D' || name[after1] == 'D' )) || ((name[before] == 'K' || name[before1] == 'K') && (name[after] == 'E' || name [after1] == 'E')) || ((name[before] == 'R' || name[before1] == 'R') && (name[after] == 'E' || name [after1] == 'E')) || ((name[before] == 'R' || name[before1] == 'R') && (name[after] == 'D' || name [after1] == 'D')) || ((name[before] == 'D' || name[before1]=='D') && (name[after] == 'R' || name[after1] == 'R')) || ((name[before] == 'D' || name[before1]=='D') && (name[after] == 'K' || name[after1]=='K')) || ((name[before] == 'E'||name[before1]=='E') && (name[after] == 'R' || name[after1]=='R')) || ((name[before] == 'E' || name[before1]=='E') && (name[after] == 'K' || name[after1]=='K'))){
        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            if(i!=1 && i!=13 && i!=2 && i!=4){
                newpopulation.push_back(name);
            }
        }
    }else if((name[after] == 'K' || name[after] == 'R' || name[after1]=='K' || name[after1]=='R') && (name[before] != 'D' && name[before] != 'E' && name[before1] != 'D' && name[before1] != 'E' ) ){
        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            if(i!=1 && i!=13){
                newpopulation.push_back(name);
            }
        }

    }else if((name[before] == 'K' || name[before] == 'R' || name[before1]=='K' || name[before1]=='R') && (name[after] != 'D' && name[after] != 'E' && name[after1] != 'D' && name[after1] != 'E') ){
        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            if(i!=1 && i!=13){
                newpopulation.push_back(name);
            }
        }

    }else if((name[before] == 'D' || name[before] == 'E' || name[before1]=='D' || name[before1]=='E') && (name[after] != 'R' && name[after] != 'K' && name[after1] != 'R' && name[after1] != 'K') ){
        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            if(i!=2 && i!=4){
                newpopulation.push_back(name);
            }
        }
    }else if((name[after] == 'D' || name[after] == 'E' || name[after1]=='D' || name[after1]=='E') && (name[before] != 'R' && name[before] != 'K' && name[before1] != 'R' && name[before1] != 'K' ) ){
        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            if(i!=2 && i!=4){
                newpopulation.push_back(name);
            }
        }
    }else if(((name[before]=='K' || name[before]=='R') && (name[before1]=='D' || name[before1]=='E')) || ((name[after]=='K' || name[after]=='R') && (name[after1]=='D' || name[after1]=='E')) || ((name[before]=='D' || name[before]=='E') && (name[before1]=='K' || name[before1]=='R')) || ((name[after]=='D' || name[after]=='E') && (name[after1]=='K' || name[after1]=='R'))){
        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            if(i!=1 && i!=13 && i!=2 && i!=4){
                newpopulation.push_back(name);
            }
        }
    }else if(((name[before]=='K' || name[before]=='R') && (name[after1]=='D' || name[after1]=='E')) || ((name[after]=='K' || name[after]=='R') && (name[before1]=='D' || name[before1]=='E')) || ((name[before]=='D' || name[before]=='E') && (name[after1]=='K' || name[after1]=='R')) || ((name[after]=='D' || name[after]=='E') && (name[before1]=='K' || name[before1]=='R'))){
        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            if(i!=1 && i!=13 && i!=2 && i!=4){
                newpopulation.push_back(name);
            }
        }
    }else{

        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            newpopulation.push_back(name);

        }
    }

}

void mutation_linear(vector<string>& newpopulation,  string name, Peptide resids){

    srand (time(NULL));
    int mut_point=rand() %20;
    int before=mut_point-1;
    int after=mut_point+1;
    if(name[before] == 'K' || name[before] == 'R' || name[after] == 'K' || name[after] == 'R'){
        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            if(i!=1 && i!=13){
                newpopulation.push_back(name);
            }
        }
    }else if(name[before] == 'D' || name[before] == 'E' || name[after] == 'D' || name[after] == 'E'){
        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            if(i!=2 && i!=4){
                newpopulation.push_back(name);
            }
        }
    }else{

        for (int i=0;i<resids.size();i++) {
            name[mut_point]=resids.name[i];
            newpopulation.push_back(name);
        }
    }
}

void crossover(Peptide pep1, Peptide pep2, vector<string>& population){
    srand (time(NULL));
    Peptide new1, new2;
    int cross_point1=rand() %20;
    int cross_point2=rand() %20;

    new1.name =pep1.name;
    new2.name=pep2.name;

    if(cross_point2 < cross_point1){
        int temp0=cross_point1;
        cross_point1=cross_point2;
        cross_point2=temp0;
    }
//cout<<new1.name<<" "<<new2.name<<endl;
    for (int i = 0; i < cross_point1; i++) {
        char temp = new1.name[i];
        new1.name[i] = new2.name[i];
        new2.name[i]= temp;
    }

    for (int i = cross_point2; i < 20; i++) {
        char temp = new1.name[i];
        new1.name[i] = new2.name[i];
        new2.name[i]= temp;
    }
    //cout<<cross_point1<<" "<<cross_point2<<endl;
    //cout<<new1.name<<" "<<new2.name<<endl;
    population.push_back(new1.name);
    population.push_back(new2.name);

}



void avg_charged(Energy_calculator calc, string input0, string inputCH, vector<vector<double>>& newEh, int type){
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

            charged_res.seq[j].init(calc.Charged[j], type);

        }
        charged_res.seq[0].E_charged_init(charged_res.seq[1], newEh, type);
}

void load_charged(Energy_calculator calc, Peptide& resids, map<char, int> res_id, int type){

    if (type==1){

        vector<vector<double>> newEh;
        avg_charged(calc, "K0h", "K+h", newEh, type);
        resids.seq[res_id['K']].init(newEh, type);

        newEh.clear();

        avg_charged(calc, "R0h", "R+h", newEh, type);
        resids.seq[res_id['R']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "D0h", "D-h", newEh, type);
        resids.seq[res_id['D']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "E0h", "E-h", newEh, type);
        resids.seq[res_id['E']].init(newEh, type);
        newEh.clear();

        /*avg_charged(calc, "H0h", "H+h", newEh);
    resids.seq[res_id['H']].init(newEh);
    newEh.clear();*/
    }
    if (type==2){
        vector<vector<double>> newEh;
        avg_charged(calc, "K0b", "K+b", newEh, type);
        resids.seq[res_id['K']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "R0b", "R+b", newEh, type);
        resids.seq[res_id['R']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "D0b", "D-b", newEh, type);
        resids.seq[res_id['D']].init(newEh, type);
        newEh.clear();

        avg_charged(calc, "E0b", "E-b", newEh, type);
        resids.seq[res_id['E']].init(newEh, type);
        newEh.clear();

        /*avg_charged(calc, "H0b", "H+b", newEh);
    resids.seq[res_id['H']].init(newEh);
    newEh.clear();*/
    }
}

void order(vector<double>& arr, vector<double>& dh,vector<double>& db, vector<string>& names)
{
    int n = arr.size();
    //sorting - ASCENDING ORDER
    for(int i=0;i<n;i++)
    {
        for(int j=i+1;j<n;j++)
        {
            if(arr[i]>arr[j])
            {
                double temp  =arr[i];
                arr[i]=arr[j];
                arr[j]=temp;
                double tempb  =db[i];
                db[i]=db[j];
                db[j]=tempb;
                double temph  =dh[i];
                dh[i]=dh[j];
                dh[j]=temph;
                string tempn  =names[i];
                names[i]=names[j];
                names[j]=tempn;
            }
        }
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
    //
    // Load energy files into resids.seq.E
    //
    for (int j=0;j<resids.name.size();j++) {
        calc.load_EfilesH(resids.name[j], j);
        calc.load_EfilesB(resids.name[j], j);
        resids.seq[j].init(calc.E_H[res_id[resids.name[j]]], 1);
        resids.seq[j].init(calc.E_B[res_id[resids.name[j]]], 2);
    }

    //
    // Load charged residues profile --> average between charged and neutral forms
    //
    load_charged(calc, resids, res_id, 1);
    load_charged(calc, resids, res_id, 2);

    //
    // Load COM of each amino acid
    //
    resids.com_aa("all.gro", res_id);

    pep.load_sequences("input");
    double maxdG=-9999999;
    double max2dG=-9999999;
    double max3dG=-9999999;
    double max4dG=-9999999;

    Peptide fittest_pep, second_fittest, third_fittest, fourth_fittest;
    int start=pep.population.size();

    /*for (int a=0;a<pep.population.size();a++) {
        pep.seq.clear();
        pep.name=pep.population[a];

        pep.seq.resize(20);
        pep.initial_pos("init_pos20", resids, res_id);
        calc.print_Emap(pep, resids, res_id, pep.name);
        calc.print_Gplot(pep, resids, res_id, pep.name);

    }
    exit(0);*/

    for (int k=0;k<start;++k) {
        mutation(pep.population, pep.population[k], resids);
        //mutation_linear(pep.population, pep.population[k], resids);
    }

    for (int i=0;i<pep.population.size();i++) {
        pep.seq.clear();
        pep.name=pep.population[i];

        //
        // Load amino acid coordinates --> alpha helix of 20 aa
        //
        pep.seq.resize(20);
        pep.initial_pos("init_pos20", resids, res_id);
        //pep.linear_pos();

        //
        // Calculate Emin for human and bacterial membrane --> stored in energy_h, energy_b
        //

        calc.get_energy_total(pep, resids, res_id);

        cout<<pep.name<<" "<< pep.energy_h<<" "<<pep.energy_b<<endl;

        //Selectivity criteria with both membranes --> maximize variable has to be max

        double maximize=pep.energy_h-pep.energy_b;

        double score=maximize;
        if(pep.depth_b>1 && pep.depth_b<3 && pep.depth_h>1 && pep.depth_h<3){
            score=score*3;
        }
        if((pep.depth_b>1 && pep.depth_b<3) || (pep.depth_h>1 && pep.depth_h<3)){
            score=score*2;
        }
        if(pep.depth_b>1 || pep.depth_h>1 ){
            score=score*2;
        }
        if (score>maxdG){
            maxdG=score;
            fittest_pep=pep;
        }

        if(score>max2dG && score<maxdG){
            max2dG=score;
            second_fittest=pep;
        }

        if(score>max3dG && score<max2dG){
            max3dG=score;
            third_fittest=pep;
        }

        if(score>max4dG && score<max2dG){
            max4dG=score;
            fourth_fittest=pep;
        }

    }


    string peptide =fittest_pep.name;

    //
    // Peptide evolution --> crossover + mutation
    //

    int rep=0;
    while(rep<10000000) {  // Mutate sequences until 100 loops happen without finding a better fit

        Peptide newpep;

            //crossover(fittest_pep,second_fittest, newpep.population);
            //crossover(third_fittest,fourth_fittest, newpep.population);
        newpep.population.push_back(fittest_pep.name);
        newpep.population.push_back(second_fittest.name);
        newpep.population.push_back(third_fittest.name);
        newpep.population.push_back(fourth_fittest.name);
            mutation(newpep.population, newpep.population[0], resids);
            mutation(newpep.population, newpep.population[1], resids);
            mutation(newpep.population, newpep.population[2], resids);
            mutation(newpep.population, newpep.population[3], resids);
            /*mutation_linear(newpep.population, newpep.population[0], resids);
            mutation_linear(newpep.population, newpep.population[1], resids);
            mutation_linear(newpep.population, newpep.population[2], resids);
            mutation_linear(newpep.population, newpep.population[3], resids);*/


        for (int i=0;i<newpep.population.size();i++) {
            //cout<<newpep.population[i]<<endl;
            newpep.seq.clear();
            newpep.name=newpep.population[i];
            newpep.seq.resize(20);
            newpep.initial_pos("init_pos20", resids, res_id);
            //newpep.linear_pos();

            calc.get_energy_total(newpep, resids, res_id);

            // selectivity criteria

                double maximize=newpep.energy_h-newpep.energy_b;

                double score=maximize;
                if(newpep.depth_b>1 && newpep.depth_b<3 && newpep.depth_h>1 && newpep.depth_h<3){
                    score=score*3;
                }
                if((newpep.depth_b>1 && newpep.depth_b<3) || (newpep.depth_h>1 && newpep.depth_h<3)){
                    score=score*1.5;
                }
                if(newpep.depth_b>1 || newpep.depth_h>1 ){
                    score=score*2;
                }
                if (score>maxdG){
                    maxdG=score;
                    fittest_pep=newpep;
                }

                if(score>max2dG && score<maxdG){
                    max2dG=score;
                    second_fittest=newpep;
                }

                if(score>max3dG && score<max2dG){
                    max3dG=score;
                    third_fittest=newpep;
                }

                if(score>max4dG && score<max2dG){
                    max4dG=score;
                    fourth_fittest=newpep;
                }

            //cout<<newpep.name<<" "<<newpep.energy_h<<" "<<newpep.energy_b<<" " <<newpep.energy_h-newpep.energy_b<<endl;

        }

        string peptide0=fittest_pep.name;

        // Sequence comparison with the last round
        if(peptide0==peptide){
            rep++;
        }else{
            rep=0;
        }

        cout<<fittest_pep.name<<" "<<rep<<" "<<fittest_pep.energy_h<<" "<<fittest_pep.energy_b<<" "<<fittest_pep.energy_h-fittest_pep.energy_b<<" "<<fittest_pep.depth_h<<" "<<fittest_pep.depth_b<<endl;

        peptide=peptide0;
        peptide0.clear();

    }

    //Print final peptide energy map and G profile
    //calc.print_Emap(fittest_pep, resids, res_id);
    //calc.print_Gplot(fittest_pep, resids, res_id);


    return 0;
}
