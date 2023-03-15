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
#include <algorithm>

#include "peptide.h"
#include "energy_calculator.h"

// 1--> human
// 2--> bacteria

using namespace std;

const double radToDeg = 57.2958;
const double degToRad = 1.0/57.2958;


void mutation_23SB(vector<string>& newpopulation,  string name, Peptide resids){

    srand (time(NULL));
    int mut_point=rand() %20;
    int before=mut_point-4;
    int before1=mut_point-3;
    int after=mut_point+4;
    int after1=mut_point+3;

    if(name[before]=='K'||name[before]=='R'||name[before1]=='K'||name[before1]=='R'){
        if(name[after]=='D'||name[after]=='E'||name[after1]=='D'||name[after1]=='E'||name[before]=='E'||name[before]=='D'||name[before1]=='D'||name[before1]=='E'){
            for (int i=0;i<resids.size()-5;i++) {
                name[mut_point]=resids.name[i];
                if(i!=1 && i!=13 && i!=2 && i!=4){
                    newpopulation.push_back(name);
                }
            }
        }else{
                for (int i=0;i<resids.size()-5;i++) {
                    name[mut_point]=resids.name[i];
                    if(i!=1 && i!=13){
                        newpopulation.push_back(name);
                    }
                }
            }
    }else if(name[before]=='E'||name[before]=='D'||name[before1]=='E'||name[before1]=='D'){
        if(name[after]=='R'||name[after]=='K'||name[after1]=='R'||name[after1]=='K'||name[before]=='R'||name[before]=='K'||name[before1]=='R'||name[before1]=='K'){
            for (int i=0;i<resids.size()-5;i++) {
                name[mut_point]=resids.name[i];
                if(i!=1 && i!=13 && i!=2 && i!=4){
                    newpopulation.push_back(name);
                }
            }
        }else{
                for (int i=0;i<resids.size()-5;i++) {
                    name[mut_point]=resids.name[i];
                    if(i!=2 && i!=4){
                        newpopulation.push_back(name);
                    }
                }
            }
    }else if(name[after]=='K'||name[after]=='R'||name[after1]=='K'||name[after1]=='R'){
        if(name[after]=='D'||name[after]=='E'||name[after1]=='D'||name[after1]=='E'||name[before]=='E'||name[before]=='D'||name[before1]=='D'||name[before1]=='E'){
            for (int i=0;i<resids.size()-5;i++) {
                name[mut_point]=resids.name[i];
                if(i!=1 && i!=13 && i!=2 && i!=4){
                    newpopulation.push_back(name);
                }
            }
        }else{
                for (int i=0;i<resids.size()-5;i++) {
                    name[mut_point]=resids.name[i];
                    if(i!=1 && i!=13){
                        newpopulation.push_back(name);
                    }
                }
            }
    }else if(name[after]=='E'||name[after]=='D'||name[after1]=='E'||name[after1]=='D'){
        if(name[after]=='R'||name[after]=='K'||name[after1]=='R'||name[after1]=='K'||name[before]=='R'||name[before]=='K'||name[before1]=='R'||name[before1]=='K'){
            for (int i=0;i<resids.size()-5;i++) {
                name[mut_point]=resids.name[i];
                if(i!=1 && i!=13 && i!=2 && i!=4){
                    newpopulation.push_back(name);
                }
            }
        }else{
                for (int i=0;i<resids.size()-5;i++) {
                    name[mut_point]=resids.name[i];
                    if(i!=2 && i!=4){
                        newpopulation.push_back(name);
                    }
                }
            }
    }else{
        for (int i=0;i<resids.size()-5;i++) {
            name[mut_point]=resids.name[i];
            newpopulation.push_back(name);

        }
    }

}

void mutation_3SB_avoid_overlap(vector<string>& newpopulation,  string name, Peptide resids){

    srand (time(NULL));
    int mut_point=rand() %20;
    int before=mut_point-3;
    int before1=mut_point-6;
    int after=mut_point+3;
    int after1=mut_point+6;

    if((name[before]=='K'||name[before]=='R') && (name[before1]=='E'||name[before1]=='D')){
            for (int i=0;i<resids.size()-11;i++) {
                name[mut_point]=resids.name[i];
                if(i!=1 && i!=13 ){
                    newpopulation.push_back(name);
                }
            }

    }else if((name[before]=='E'||name[before]=='D') && (name[before1]=='K'||name[before1]=='R')){
            for (int i=0;i<resids.size()-11;i++) {
                name[mut_point]=resids.name[i];
                if( i!=2 && i!=4){
                    newpopulation.push_back(name);
                }
            }
    }else if((name[after]=='K'||name[after]=='R') && (name[after1]=='E'||name[after1]=='D')){
            for (int i=0;i<resids.size()-11;i++) {
                name[mut_point]=resids.name[i];
                if(i!=1 && i!=13 ){
                    newpopulation.push_back(name);
                }
            }
    }else if((name[after]=='D'||name[after]=='E') && (name[after1]=='K'||name[after1]=='R')){
            for (int i=0;i<resids.size()-11;i++) {
                name[mut_point]=resids.name[i];
                if(i!=2 && i!=4){
                    newpopulation.push_back(name);
                }
            }
    }else if((name[after]=='D'||name[after]=='E') && (name[before]=='D'||name[before]=='E')){
        for (int i=0;i<resids.size()-11;i++) {
            name[mut_point]=resids.name[i];
            if(i!=2 && i!=4){
                newpopulation.push_back(name);
            }
        }
    }else if((name[after]=='K'||name[after]=='R') && (name[before]=='K'||name[before]=='R')){
        for (int i=0;i<resids.size()-11;i++) {
            name[mut_point]=resids.name[i];
            if(i!=1 && i!=13){
                newpopulation.push_back(name);
            }
        }
    }else if((name[after]=='W'|| name[after]=='Y' || name[after]=='F') && (name[before]=='W'||name[before]=='Y' || name[before]=='F')){
        for (int i=0;i<resids.size()-11;i++) {
            name[mut_point]=resids.name[i];
            if(i!=3 && i!=11 && i!=16){
                newpopulation.push_back(name);
            }
        }
    }else if((name[after]=='W'|| name[after]=='Y' || name[after]=='F') && (name[after1]=='W'||name[after1]=='Y' || name[after1]=='F')){
        for (int i=0;i<resids.size()-11;i++) {
            name[mut_point]=resids.name[i];
            if(i!=3 && i!=11 && i!=16){
                newpopulation.push_back(name);
            }
        }
    }else if((name[before]=='W'|| name[before]=='Y' || name[before]=='F') && (name[before1]=='W'||name[before1]=='Y' || name[before1]=='F')){
        for (int i=0;i<resids.size()-11;i++) {
            name[mut_point]=resids.name[i];
            if(i!=3 && i!=11 && i!=16){
                newpopulation.push_back(name);
            }
        }
    }else{
        for (int i=0;i<resids.size()-11;i++) {
            name[mut_point]=resids.name[i];
            newpopulation.push_back(name);

        }
    }

}

bool salt_bridges(Peptide pep){
    int SB=0;
    for (int i=0;i<pep.seq.size()-3;++i){
        if((pep.name[i] == 'K' || pep.name[i] == 'R') && (pep.name[i+3] == 'D' || pep.name[i+3] == 'E')){
            SB=SB+1;
        }else if((pep.name[i] == 'D' || pep.name[i] == 'E') && (pep.name[i+3] == 'K' || pep.name[i+3] == 'R')){
            SB=SB+1;
        }else{
            SB=SB+0;
        }
    }
    if(SB>0){
        return true;
    }else{
        return false;
    }
}

bool aromatics(Peptide pep){
    int AR=0;
    for (int i=0;i<pep.seq.size()-3;++i){
        if(pep.name[i] == 'W' && pep.name[i+3] == 'W'){
            AR=AR+1;
        }if(pep.name[i] == 'Y' && pep.name[i+3] == 'Y'){
            AR=AR+1;
        }if(pep.name[i] == 'F' && pep.name[i+3] == 'F'){
            AR=AR+1;
        }if((pep.name[i] == 'W' && pep.name[i+3] == 'F') || (pep.name[i] == 'F' && pep.name[i+3] == 'W')){
            AR=AR+1;
        }if((pep.name[i] == 'W' && pep.name[i+3] == 'Y') || (pep.name[i] == 'Y' && pep.name[i+3] == 'W')){
            AR=AR+1;
        }else if((pep.name[i] == 'Y' && pep.name[i] == 'F') || (pep.name[i+3] == 'F' && pep.name[i+3] == 'Y')){
            AR=AR+1;
        }else{
            AR=AR+0;
        }
    }
    if(AR>0){
        return true;
    }else{
        return false;
    }
}

bool find_overlapping_SB(string name){
    vector<int>index1;
    vector<int>index2;

    for(int i=0;i<name.size()-3;++i){
        if((name[i] == 'K' && name[i+3] == 'D') || (name[i] == 'D' && name[i+3] == 'K')){
            index1.push_back(i);
            index2.push_back(i+3);
        }
        if((name[i] == 'K' &&  name[i+3] == 'E') || (name[i] == 'E' && name[i+3] == 'K')){
            index1.push_back(i);
            index2.push_back(i+3);
        }
        if((name[i] == 'R' && name[i+3] == 'D' ) || (name[i] == 'D' && name[i+3] == 'R')){
            index1.push_back(i);
            index2.push_back(i+3);
        }
        if(( name[i] == 'R' &&  name[i+3] == 'E') || (name[i] == 'E' && name[i+3] == 'R')){
            index1.push_back(i);
            index2.push_back(i+3);
        }
    }

    //Create batch of possible SB interactions and change the COM of the SB
    vector<int>rmove_idx;
    for(int s=0;s<index1.size();++s){
        for(int t=0;t<index1.size();++t){
            if(index1[s]==index2[t]){
                rmove_idx.push_back(s);
            }
        }
        //cout<<index1[s]<<" "<<index2[s]<<endl;
    }

    if(rmove_idx.size()>0){
        index1.clear();
        index2.clear();
        return true;
    }else{
        index1.clear();
        index2.clear();
        return false;
    }
}

void detect_interactions(Peptide pep, std::vector<int>&index1, std::vector<int>&index2, vector<char>&SB_type){
    for(int i=0;i<pep.name.size()-3;++i){
        if((pep.name[i] == 'K' && pep.name[i+3] == 'D') || (pep.name[i] == 'D' && pep.name[i+3] == 'K')){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('1');
        }
        if((pep.name[i] == 'K' &&  pep.name[i+3] == 'E') || (pep.name[i] == 'E' && pep.name[i+3] == 'K')){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('2');
        }
        if((pep.name[i] == 'R' && pep.name[i+3] == 'D' ) || (pep.name[i] == 'D' && pep.name[i+3] == 'R')){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('3');
        }
        if(( pep.name[i] == 'R' &&  pep.name[i+3] == 'E') || (pep.name[i] == 'E' && pep.name[i+3] == 'R')){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('4');
        }
        if(pep.name[i] == 'W' && pep.name[i+3] == 'W'){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('5');
        }
        if(pep.name[i] == 'F' && pep.name[i+3] == 'F'){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('6');
        }
        if(pep.name[i] == 'Y' && pep.name[i+3] == 'Y'){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('7');
        }
        if((pep.name[i] == 'W' &&  pep.name[i+3] == 'F') || (pep.name[i] == 'F' && pep.name[i+3] == 'W')){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('8');
        }
        if((pep.name[i] == 'W' && pep.name[i+3] == 'Y' ) || (pep.name[i] == 'Y' && pep.name[i+3] == 'W')){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('9');
        }
        if(( pep.name[i] == 'Y' &&  pep.name[i+3] == 'F') || (pep.name[i] == 'F' && pep.name[i+3] == 'Y')){
            index1.push_back(i);
            index2.push_back(i+3);
            SB_type.push_back('0');
        }
    }
}


void find_index_overlap(vector<int>&rmove_idx,vector<int>&rmove_idx2, std::vector<int>index1, std::vector<int>index2){
    for(int s=0;s<index1.size();++s){
        for(int t=0;t<index1.size();++t){
            if(index1[s]==index2[t]){
                rmove_idx.push_back(s);
                rmove_idx2.push_back(t);
            }
        }
    }
}

void make_subtitutions(Peptide prev_pep,  std::vector<int>index1, std::vector<int>index2, vector<char>SB_type,  Peptide& new_pep){

    int count=0;
    int p=0;
    new_pep.seq.resize(prev_pep.seq.size()-index1.size());

    for(int k=0;k<prev_pep.name.size();++k){

        if (std::find(index1.begin(), index1.end(), k) != index1.end()) {
            new_pep.name.push_back(SB_type[count]);
            new_pep.seq[p].pos = (prev_pep.seq[k].pos + prev_pep.seq[k+3].pos)/2;
            count+=1;
            p+=1;

        }
        else if((std::find(index2.begin(), index2.end(), k) != index2.end()) || (std::find(index1.begin(), index1.end(), k) != index1.end())){
            new_pep.name.push_back('X');

        }else{
            new_pep.name.push_back(prev_pep.name[k]);
            new_pep.seq[p].pos=prev_pep.seq[k].pos;
            p+=1;
        }

    }
}

void rename(Peptide& pep1){
    string newname;
    for (int x=0;x<pep1.name.size();++x){
        if (pep1.name[x] != 'X'){
            newname=newname+pep1.name[x];

        }
    }
    pep1.name.clear();
    pep1.name=newname;
}

void population_SB(  Peptide pep, vector<Peptide>&SB_all){
    vector<char>SB_type;
    std::vector<int>index1;
    std::vector<int>index2;

    //Detect salt-bridges, store their position and type
    detect_interactions(pep,index1,index2,SB_type);

    //Detect if salt-bridges overlap; if there is overlap, rmove_idx will be bigger than 0
    vector<int>rmove_idx;
    vector<int>rmove_idx2;
    find_index_overlap(rmove_idx,rmove_idx2,index1,index2);

    //Create batch of possible SB interactions and change the COM of the SB
    if( rmove_idx.size() == 0 ){

        Peptide pep1;

        make_subtitutions(pep, index1, index2, SB_type, pep1);
        rename(pep1);

        SB_all.push_back(pep1);

    }else{

        cout<<"overlapping interactions: "<< pep.name<<endl;
        exit(1);
        vector<int>new_idx1;
        vector<int>new_idx2;
        for(int i=0;i<rmove_idx.size();++i){
            cout<<rmove_idx[i]<<" "<<rmove_idx2[i]<<endl;

        }

        /*for(int j=0;j<index1.size();++j){
                //cout<<rmove_idx[i] <<" "<<j<<endl;
                if ( ! (std::find(rmove_idx.begin(), rmove_idx.end(), j) != rmove_idx.end())){
                    new_idx1.push_back(index1[j]);
                    new_idx2.push_back(index2[j]);
                }
                //cout<<new_idx1.back()<<endl;
                if ( ! (std::find(rmove_idx2.begin(), rmove_idx2.end(), j) != rmove_idx2.end()) ){
                   /* if (new_idx1.size() == 0){
                    new_idx1.push_back(index1[j]);
                    new_idx2.push_back(index2[j]);
                    }else{
                        //if(new_idx1.back() != index1[j]){
                            new_idx1.push_back(index1[j]);
                            new_idx2.push_back(index2[j]);
                        //}
                    //}
                }

        }*/



            int count=0;
            for(int a=0;a<new_idx1.size();++a){
                //cout<<i<<endl;
                Peptide new_pep;

                int p=0;
                new_pep.seq.resize(19);
                for(int k=0;k<pep.name.size();++k){

                    if (new_idx1[a] == k){
                        //if (std::find(new_idx1.begin(), new_idx1.end(), k) != new_idx1.end()) {
                        new_pep.name.push_back(SB_type[count]);
                        new_pep.seq[p].pos = (pep.seq[k].pos + pep.seq[k+3].pos)/2;
                        count+=1;
                        p+=1;

                    }else if (new_idx2[a] == k){
                        new_pep.name.push_back('X');
                    }else{
                        new_pep.name.push_back(pep.name[k]);
                        new_pep.seq[p].pos=pep.seq[k].pos;
                        p+=1;
                    }

                }

                cout<<"new_pep:" <<new_pep.name<<endl;
                SB_all.push_back(new_pep);
            }
    }

}

void mutation_linear(vector<string>& newpopulation,  string name, Peptide resids){

    srand (time(NULL));
    int mut_point=rand() %20;
    int before=mut_point-1;
    int after=mut_point+1;
    if(name[before] == 'K' || name[before] == 'R' || name[after] == 'K' || name[after] == 'R'){
        for (int i=0;i<resids.size()-5;i++) {
            name[mut_point]=resids.name[i];
            if(i!=1 && i!=13){
                newpopulation.push_back(name);
            }
        }
    }else if(name[before] == 'D' || name[before] == 'E' || name[after] == 'D' || name[after] == 'E'){
        for (int i=0;i<resids.size()-5;i++) {
            name[mut_point]=resids.name[i];
            if(i!=2 && i!=4){
                newpopulation.push_back(name);
            }
        }
    }else{

        for (int i=0;i<resids.size()-5;i++) {
            name[mut_point]=resids.name[i];
            newpopulation.push_back(name);
        }
    }
}

void mutation_NR(vector<string>& newpopulation,  string name, Peptide resids){
    srand (time(NULL));
    int mut_point=rand() %20;
        for (int i=0;i<resids.size()-5;i++) {
            name[mut_point]=resids.name[i];
            newpopulation.push_back(name);
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

    Peptide resids;
    resids.name="AEKWRCVLIMQFNDSTY1234567890O";
    resids.seq.resize(resids.name.size());
    map<char, int> res_id{ { 'A', 0 }, { 'E', 1 }, { 'K', 2 },  { 'W', 3 }, { 'R', 4 }, { 'C', 5 }, { 'V', 6 }, { 'L', 7 }, { 'I', 8 },  { 'M', 9 }, { 'Q', 10 }, { 'F', 11 },  { 'N', 12 }, {'D', 13}, {'S', 14}, {'T', 15}, {'Y', 16},{'1',17},{'2',18},{'3',19},{'4',20},{'5',21},{'6',22},{'7',23},{'8',24},{'9',25},{'0',26},{'O',27}};

    calc.init_grid(resids.name.size(), 200, 2 );
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
        //pep.seq.clear();
        pep.name=pep.population[a];

        pep.seq.resize(20);
        pep.initial_pos("init_pos20", resids, res_id);
        //pep.linear_pos();
        if(salt_bridges(pep) || aromatics(pep)){

            vector<Peptide>SB_all;

            population_SB(pep, SB_all);
            calc.print_Gplot_SB(SB_all,resids,res_id, pep);
            calc.print_Gplot(pep, resids, res_id, pep.name);
           // calc.print_Emap(SB_all[0], resids, res_id, SB_all[0].name);
        }else{
            calc.print_Emap(pep, resids, res_id, pep.name);
            calc.print_Gplot(pep, resids, res_id, pep.name);
        }

    }
    exit(0);*/

    for (int k=0;k<start;++k) {
        //mutation_NR(pep.population, pep.population[k], resids);
        //if (find_overlapping_SB(pep.population[k])){
            mutation_3SB_avoid_overlap(pep.population, pep.population[k], resids);
        /*}else {
                mutation_NR(pep.population, pep.population[k], resids);
        }*/
        //mutation_23SB(pep.population, pep.population[k], resids);
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

        if(salt_bridges(pep) || aromatics(pep)){

            vector<Peptide>SB_all;

            population_SB( pep, SB_all);
            //calculate energy of each of the peptides and do the avg along the z depth.

            calc.get_energy_SB(SB_all,resids,res_id, pep);
            cout<<"SB "<<pep.name<<" "<<pep.energy_h<<" "<<pep.energy_b<<" "<<pep.energy_h-pep.energy_b<<" "<<pep.depth_h<<" "<<pep.depth_b<<endl;

        }else{

            calc.get_energy_total(pep, resids, res_id);

            cout<<pep.name<<" "<< pep.energy_h<<" "<<pep.energy_b<<" "<<pep.energy_h-pep.energy_b<<" "<<pep.depth_h<<" "<<pep.depth_b<<endl;
        }

        //Selectivity criteria with both membranes --> maximize variable has to be max

        double maximize=pep.energy_h-pep.energy_b;

        double score=maximize;
        /*if(pep.depth_b>1 && pep.depth_b<3 && pep.depth_h>1 && pep.depth_h<3){
            score=score*4;
        }
        if((pep.depth_b>1 && pep.depth_b<3) || (pep.depth_h>1 && pep.depth_h<3)){
            score=score*3;
        }
        if(pep.depth_b>1 || pep.depth_h>1 ){
            score=score*3;
        }*/
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

    ofstream file1;
    file1.open("sequences+50");
    ofstream file2;
    file2.open("sequences-39h");

    int rep=0;
    while(rep<10000) {  // Mutate sequences until 100 loops happen without finding a better fit

        Peptide newpep;

            //crossover(fittest_pep,second_fittest, newpep.population);
            //crossover(third_fittest,fourth_fittest, newpep.population);
        newpep.population.push_back(fittest_pep.name);
        newpep.population.push_back(second_fittest.name);
        newpep.population.push_back(third_fittest.name);
        newpep.population.push_back(fourth_fittest.name);
        int start2=newpep.population.size();

        for (int k=0;k<start2;++k){
            //if (find_overlapping_SB(newpep.population[k])){
                mutation_3SB_avoid_overlap(newpep.population, newpep.population[k], resids);
            /*}else {
                mutation_NR(newpep.population, newpep.population[k], resids);
            }*/
        }


        for (int i=0;i<newpep.population.size();i++) {
            newpep.seq.clear();
            newpep.name=newpep.population[i];
            newpep.seq.resize(20);
            newpep.initial_pos("init_pos20", resids, res_id);
            //newpep.linear_pos();

            if(salt_bridges(newpep) || aromatics(newpep)){

                vector<Peptide>SB_all;

                population_SB(  newpep, SB_all);
                calc.get_energy_SB(SB_all,resids,res_id, newpep);

            }else{

                calc.get_energy_total(newpep, resids, res_id);

            }

            // selectivity criteria

                double maximize=newpep.energy_h-newpep.energy_b;

                double score=maximize;
                /*if(pep.depth_b>1 && pep.depth_b<3 && pep.depth_h>1 && pep.depth_h<3){
                    score=score*4;
                }
                if((pep.depth_b>1 && pep.depth_b<3) || (pep.depth_h>1 && pep.depth_h<3)){
                    score=score*3;
                }
                if(pep.depth_b>1 || pep.depth_h>1 ){
                    score=score*3;
                }*/

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

        if(rep==0 && fittest_pep.energy_h-fittest_pep.energy_b > 50){
            file1<<fittest_pep.name<<" "<<rep<<" "<<fittest_pep.energy_h<<" "<<fittest_pep.energy_b<<" "<<fittest_pep.energy_h-fittest_pep.energy_b<<" "<<fittest_pep.depth_h<<" "<<fittest_pep.depth_b<<endl;
        }
        if(rep==0 && fittest_pep.depth_h>3.3 && fittest_pep.depth_b < 3){
            file2<<fittest_pep.name<<" "<<rep<<" "<<fittest_pep.energy_h<<" "<<fittest_pep.energy_b<<" "<<fittest_pep.energy_h-fittest_pep.energy_b<<" "<<fittest_pep.depth_h<<" "<<fittest_pep.depth_b<<endl;
        }

        peptide=peptide0;
        peptide0.clear();

    }
    file1.close();
    file2.close();

    //Print final peptide energy map and G profile
    //calc.print_Emap(fittest_pep, resids, res_id);
    //calc.print_Gplot(fittest_pep, resids, res_id);


    return 0;
}
