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
    std::vector<int> list{0,2,3,6,7,9,10,13,14,16,17};
    int index = rand() % list.size(); // pick a random index
    int mut_point = list[index];

    //int mut_point=rand() %20;
    int before=mut_point-3;
    int before1=mut_point-6;
    int after=mut_point+3;
    int after1=mut_point+6;
    //int i=rand() %16;
    //newpopulation.push_back(name);

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

void add_mut(vector<string>& newpopulation,  string name, Peptide resids){

    //srand (time(NULL));
    int mut_point=rand() %20;
    int i=rand() %17;
    int before=mut_point-3;
    int before1=mut_point-6;
    int after=mut_point+3;
    int after1=mut_point+6;

    if((name[before]=='K'||name[before]=='R') && (name[before1]=='E'||name[before1]=='D')){

                name[mut_point]=resids.name[i];
                if(i!=1 && i!=13 ){
                    newpopulation.push_back(name);
                }


    }else if((name[before]=='E'||name[before]=='D') && (name[before1]=='K'||name[before1]=='R')){

                name[mut_point]=resids.name[i];
                if( i!=2 && i!=4){
                    newpopulation.push_back(name);
                }

    }else if((name[after]=='K'||name[after]=='R') && (name[after1]=='E'||name[after1]=='D')){

                name[mut_point]=resids.name[i];
                if(i!=1 && i!=13 ){
                    newpopulation.push_back(name);
                }

    }else if((name[after]=='D'||name[after]=='E') && (name[after1]=='K'||name[after1]=='R')){
                name[mut_point]=resids.name[i];
                if(i!=2 && i!=4){
                    newpopulation.push_back(name);
                }

    }else if((name[after]=='D'||name[after]=='E') && (name[before]=='D'||name[before]=='E')){
            name[mut_point]=resids.name[i];
            if(i!=2 && i!=4){
                newpopulation.push_back(name);
            }

    }else if((name[after]=='K'||name[after]=='R') && (name[before]=='K'||name[before]=='R')){
            name[mut_point]=resids.name[i];
            if(i!=1 && i!=13){
                newpopulation.push_back(name);
            }

    }else if((name[after]=='W'|| name[after]=='Y' || name[after]=='F') && (name[before]=='W'||name[before]=='Y' || name[before]=='F')){
            name[mut_point]=resids.name[i];
            if(i!=3 && i!=11 && i!=16){
                newpopulation.push_back(name);
            }

    }else if((name[after]=='W'|| name[after]=='Y' || name[after]=='F') && (name[after1]=='W'||name[after1]=='Y' || name[after1]=='F')){
            name[mut_point]=resids.name[i];
            if(i!=3 && i!=11 && i!=16){
                newpopulation.push_back(name);
            }

    }else if((name[before]=='W'|| name[before]=='Y' || name[before]=='F') && (name[before1]=='W'||name[before1]=='Y' || name[before1]=='F')){
            name[mut_point]=resids.name[i];
            if(i!=3 && i!=11 && i!=16){
                newpopulation.push_back(name);
            }

    }else{
            name[mut_point]=resids.name[i];
            newpopulation.push_back(name);

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
        }else if((pep.name[i] == 'Y' && pep.name[i+3] == 'F') || (pep.name[i] == 'F' && pep.name[i+3] == 'Y')){
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
        if(name[i] == 'W' && name[i+3] == 'W'){
            index1.push_back(i);
            index2.push_back(i+3);
        }
        if(name[i] == 'F' && name[i+3] == 'F'){
            index1.push_back(i);
            index2.push_back(i+3);
        }
        if(name[i] == 'Y' && name[i+3] == 'Y'){
            index1.push_back(i);
            index2.push_back(i+3);
        }
        if((name[i] == 'W' &&  name[i+3] == 'F') || (name[i] == 'F' && name[i+3] == 'W')){
            index1.push_back(i);
            index2.push_back(i+3);
        }
        if((name[i] == 'W' && name[i+3] == 'Y' ) || (name[i] == 'Y' && name[i+3] == 'W')){
            index1.push_back(i);
            index2.push_back(i+3);
        }
        if(( name[i] == 'Y' &&  name[i+3] == 'F') || (name[i] == 'F' && name[i+3] == 'Y')){
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
    //count=rmove_idx.size();

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

void detect_interactions(Peptide pep, vector<int>&index1, vector<int>&index2, vector<char>&SB_type){
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


void find_index_overlap(vector<int>&rmove_idx,vector<int>&rmove_idx2, vector<int>&index1, vector<int>&index2, vector<char>&SB_type, vector<char>&SB_type_overlap){
    for(int s=0;s<index1.size();++s){
        for(int t=0;t<index1.size();++t){
            if(index1[s]==index2[t]){
                rmove_idx.push_back(s);
                rmove_idx2.push_back(t);
                SB_type_overlap.push_back(SB_type[t]);
                SB_type_overlap.push_back(SB_type[s]);

                SB_type.erase (SB_type.begin()+s);
                SB_type.erase (SB_type.begin()+t);
            }
        }
    }
}

void find_index_overlap_cross(vector<int>&rmove_idx,vector<int>&rmove_idx2, vector<int>&index1, vector<int>&index2, vector<char>&SB_type, vector<char>&SB_type_overlap){
    for(int s=0;s<index1.size();++s){
        for(int t=0;t<index1.size();++t){
            if(index1[s]==index2[t]){
                rmove_idx.push_back(s);
                rmove_idx2.push_back(t);
                SB_type_overlap.push_back(SB_type[t]);
                SB_type_overlap.push_back(SB_type[s]);

            }
        }
    }
}

void make_subtitutions(Peptide&prev_pep,  vector<int>index1, vector<int>index2, vector<char>SB_type,  Peptide& new_pep){

    int count=0;
    int p=0;
    int r=0;
  //  new_pep.seq.resize(prev_pep.seq.size()-index1.size());

    //Change residues positions

    for(int k=0;k<prev_pep.name.size();++k){

        if (std::find(index1.begin(), index1.end(), k) != index1.end()) {
            new_pep.seq[p].pos = (prev_pep.seq[r].pos + prev_pep.seq[r+3].pos)/2;
            p+=1;
            r+=1;
        }
        else if((std::find(index2.begin(), index2.end(), k) != index2.end()) || (std::find(index1.begin(), index1.end(), k) != index1.end()) ){
            r+=1;

        }else if( prev_pep.name[k]=='X'){

        }else{
            new_pep.seq[p].pos=prev_pep.seq[r].pos;
            p+=1;
            r+=1;
        }
    }

    //Change peptide name

    for(int k=0;k<prev_pep.name.size();++k){

        if (std::find(index1.begin(), index1.end(), k) != index1.end()) {
            new_pep.name.push_back(SB_type[count]);
            count+=1;
        }
        else if((std::find(index2.begin(), index2.end(), k) != index2.end()) || (std::find(index1.begin(), index1.end(), k) != index1.end())){
            new_pep.name.push_back('X');

        }else{
            new_pep.name.push_back(prev_pep.name[k]);
        }
    }
}

void make_substitutions_2(Peptide&prev_pep,  int n1, int n2, vector<char>SB_type,  Peptide& new_pep, int& count){

    int p=0;
  //  new_pep.seq.resize(prev_pep.seq.size()-1);

    //Change peptide name

    for(int k=0;k<prev_pep.name.size();++k){

        if ( k == n1) {
            new_pep.name.push_back(SB_type[count]);
            count+=1;
        }
        else if( k == n2){
            new_pep.name.push_back('X');

        }else{
            new_pep.name.push_back(prev_pep.name[k]);
        }

    }

    //Change residues positions

    for(int k=0;k<prev_pep.seq.size();++k){

        if ( k == n1) {
            new_pep.seq[p].pos = (prev_pep.seq[k].pos + prev_pep.seq[k+3].pos)/2;
            p+=1;

        }
        else if( k == n2){


        }else{
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
    vector<int>index1;
    vector<int>index2;
    vector<char>SB_type_overlap;

    //Detect salt-bridges, store their position and type
    detect_interactions(pep,index1,index2,SB_type);

    //Detect if salt-bridges overlap; if there is overlap, rmove_idx will be bigger than 0
    vector<int>rmove_idx;
    vector<int>rmove_idx2;
    find_index_overlap(rmove_idx,rmove_idx2,index1,index2, SB_type,SB_type_overlap);

    //Create batch of possible SB interactions and change the COM of the SB
    if( rmove_idx.size() == 0 ){

        Peptide pep1;

        pep1.seq.resize(pep.name.size()-index1.size());
        make_subtitutions(pep, index1, index2, SB_type, pep1);
        rename(pep1);

        SB_all.push_back(pep1);

    }else{


        //cout<<"overlapping interactions: "<< pep.name<<endl;

        Peptide pep2;
        pep2.seq.resize(pep.name.size()-1);
        int ct=0;
        make_substitutions_2(pep,index1[rmove_idx2[0]],index1[rmove_idx[0]],SB_type_overlap,pep2,ct);
        Peptide pep3;
        pep3.seq.resize(pep.name.size()-1);
        make_substitutions_2(pep,index2[rmove_idx2[0]],index2[rmove_idx[0]],SB_type_overlap,pep3,ct);

        index1.erase(index1.begin()+rmove_idx[0]);
        index1.erase(index1.begin()+rmove_idx2[0]);
        index2.erase(index2.begin()+rmove_idx[0]);
        index2.erase(index2.begin()+rmove_idx2[0]);


        Peptide pep21;
        Peptide pep31;
        pep21.seq.resize(pep2.name.size()-index1.size());
        make_subtitutions(pep2,index1,index2,SB_type,pep21);

        rename(pep21);
        SB_all.push_back(pep21);

        pep31.seq.resize(pep3.name.size()-index1.size());
        make_subtitutions(pep3,index1,index2,SB_type,pep31);
        rename(pep31);
        SB_all.push_back(pep31);

        //cout<<SB_all[0].name<<endl;
        //cout<<SB_all[1].name<<endl;

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
    std::vector<int> list{0,2,3,6,7,9,10,13,14,16,17};
    int index = rand() % list.size(); // pick a random index
    int mut_point = list[index];
    //int mut_point=rand() %20;
    int before=mut_point-3;
    int before1=mut_point-6;
    int after=mut_point+3;
    int after1=mut_point+6;
  //  int i=rand() %16;
   // newpopulation.push_back(name);

    for (int i=0;i<resids.size()-11;i++) {
        name[mut_point]=resids.name[i];

        if((name[before]=='K'||name[before]=='R') && (name[before1]=='E'||name[before1]=='D') && (name[after]=='R'||name[after]=='K')){
                if(i!=1 && i!=13 ){
                    newpopulation.push_back(name);
                }
        }else if((name[before]=='E'||name[before]=='D') && (name[before1]=='K'||name[before1]=='R') && (name[after]=='E'||name[after]=='D')){
                if( i!=2 && i!=4){
                    newpopulation.push_back(name);
                }

        }else if((name[before]=='E'||name[before]=='D') && (name[after]=='D'||name[after]=='E') && (name[after1]=='R'||name[after1]=='K')){
                if( i!=2 && i!=4){
                    newpopulation.push_back(name);
                }

        }else if((name[before]=='R'||name[before]=='K') && (name[after]=='R'||name[after]=='K') && (name[after1]=='D'||name[after1]=='E')){
                if( i!=1 && i!=13){
                    newpopulation.push_back(name);
                }

        }else if((name[before]=='W'||name[before]=='F'||name[before]=='Y') && (name[before1]=='W'||name[before1]=='F'||name[before1]=='Y') && (name[after]=='W'||name[after]=='F'||name[after]=='Y')){

                if(i!=3 && i!=11 && i!=16 ){
                    newpopulation.push_back(name);
                }
        }else if((name[after]=='W'||name[after]=='F'||name[after]=='Y') && (name[after1]=='W'||name[after1]=='F'||name[after1]=='Y') && (name[before]=='W'||name[before]=='F'||name[before]=='Y')){
                if(i!=3 && i!=11 && i!=16 ){
                    newpopulation.push_back(name);
                }
        }else{

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
    vector<char>SB_type;
    vector<int>index1;
    vector<int>index2;
    vector<char>SB_type_overlap;

    //Detect salt-bridges, store their position and type
    detect_interactions(new1,index1,index2,SB_type);

    //Detect if salt-bridges overlap; if there is overlap, rmove_idx will be bigger than 0
    vector<int>rmove_idx;
    vector<int>rmove_idx2;
    find_index_overlap_cross(rmove_idx,rmove_idx2,index1,index2, SB_type,SB_type_overlap);
    if (rmove_idx.size()<2){
        population.push_back(new1.name);
    }

    /*rmove_idx.clear();
    rmove_idx2.clear();
    SB_type.clear();
    index1.clear();
    index2.clear();
    SB_type_overlap.clear();*/


        //population.push_back(new2.name);
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

string gen_random(){
    string res="AEKWRCVLIMQFNDSTY";
    int done=0;

    while(done==0){
        Peptide random_pep;
        std::random_device rd;
        std::mt19937 g(rd());
        std::uniform_int_distribution<> dis(0, res.size() - 1);

        for (int i = 0; i < 20; ++i) {
                random_pep.name += res[dis(g)];
        }


            vector<char>SB_type;
            vector<int>index1;
            vector<int>index2;
            vector<char>SB_type_overlap;

            //Detect salt-bridges, store their position and type
            detect_interactions(random_pep,index1,index2,SB_type);

            //Detect if salt-bridges overlap; if there is overlap, rmove_idx will be bigger than 0
            vector<int>rmove_idx;
            vector<int>rmove_idx2;
            find_index_overlap_cross(rmove_idx,rmove_idx2,index1,index2, SB_type,SB_type_overlap);
            if (rmove_idx.size()<2){
                done=1;
                return random_pep.name;
            }else{
                done=0;
            }
    }

}

void random_input(int number, Peptide& curr_popu){
    while(curr_popu.population.size()<number){
         Peptide rando;
         rando.name=gen_random();
         vector<char>SB_type;
         vector<int>index1;
         vector<int>index2;
         vector<char>SB_type_overlap;

         //Detect salt-bridges, store their position and type
         detect_interactions(rando,index1,index2,SB_type);

         //Detect if salt-bridges overlap; if there is overlap, rmove_idx will be bigger than 0
         vector<int>rmove_idx;
         vector<int>rmove_idx2;
         find_index_overlap_cross(rmove_idx,rmove_idx2,index1,index2, SB_type,SB_type_overlap);
         if (rmove_idx.size()<2){
             curr_popu.population.push_back(rando.name);
         }
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

    std::fstream file( "input", std::fstream::in );
    file.seekg(0,ios::end);

    if(file.tellg()>=1){
        pep.load_sequences("input");
    }
    double maxdG=-9999999;
    double max2dG=-9999999;
    double max3dG=-9999999;
    double max4dG=-9999999;

    double mindG=9999999;
    double min2dG=9999999;
    double min3dG=9999999;
    double min4dG=9999999;


    //random_input(100, pep); //generates X number of random sequences and adds them to pep.population


    int start=pep.population.size();

    for (int a=0;a<pep.population.size();a++) {
        //pep.seq.clear();
        pep.name=pep.population[a];

        pep.seq.resize(20);
        pep.initial_pos("init_pos20", resids, res_id);
        //pep.linear_pos();

        //calc.get_energy_total(pep, resids, res_id);

        //cout<<pep.name<<" "<< pep.energy_h<<" "<<pep.energy_b<<" "<<pep.energy_h-pep.energy_b<<" "<<pep.depth_h<<" "<<pep.depth_b<<endl;

        //pep.linear_pos();
        if(salt_bridges(pep) || aromatics(pep)){

            vector<Peptide>SB_all;

            population_SB(pep, SB_all);
            calc.print_Gplot_SB(SB_all,resids,res_id, pep);
            //calc.print_Gplot(pep, resids, res_id, pep.name);
            //calc.print_Emap(pep, resids, res_id, pep.name);
            //calc.print_Emap(SB_all[0], resids, res_id, SB_all[0].name);
            calc.get_energy_SB(SB_all,resids,res_id, pep);
            cout<<pep.name<<" "<<pep.energy_h<<" "<<pep.energy_b<<" "<<pep.energy_h-pep.energy_b<<" "<<pep.depth_h<<" "<<pep.depth_b<<endl;

        }else{
            //calc.print_Emap(pep, resids, res_id, pep.name);
            calc.print_Gplot(pep, resids, res_id, pep.name);
            calc.get_energy_total(pep, resids, res_id);

            cout<<pep.name<<" "<< pep.energy_h<<" "<<pep.energy_b<<" "<<pep.energy_h-pep.energy_b<<" "<<pep.depth_h<<" "<<pep.depth_b<<endl;
        }

    }
    exit(0);

    /*for (int k=0;k<start;++k) {
        if (find_overlapping_SB(pep.population[k])){
            mutation_3SB_avoid_overlap(pep.population, pep.population[k], resids);
        }else {
                mutation_NR(pep.population, pep.population[k], resids);
        }
    }*/
    Peptide fittest_pep, second_fittest, third_fittest, fourth_fittest;
    for (int i=0;i<pep.population.size();i++) {
        //pep.seq.clear();
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
        //double minimize=pep.energy_h;
        double maximize=pep.energy_h-pep.energy_b;
        double max_z=pep.depth_h-pep.depth_b;
        double score=maximize;
        //double score=1000*max_z+0.1*maximize;
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
            fourth_fittest=third_fittest;
            third_fittest=second_fittest;
            second_fittest=fittest_pep;
            fittest_pep=pep;

        }

        if(score>max2dG && score<maxdG){
            max2dG=score;
            fourth_fittest=third_fittest;
            third_fittest=second_fittest;
            second_fittest=pep;
        }

        if(score>max3dG && score<max2dG){
            max3dG=score;
            fourth_fittest=third_fittest;
            third_fittest=pep;
        }

        if(score>max4dG && score<max3dG){
            max4dG=score;
            fourth_fittest=pep;
        }
        /*if (pep.depth_h>0.5 && pep.depth_h<1.5 && minimize<mindG){
            mindG=minimize;
            fittest_pep=pep;
        }

        if(pep.depth_h>0.5 && pep.depth_h<1.5 && minimize<min2dG && minimize>mindG){
            min2dG=minimize;
            second_fittest=pep;
        }

        if(pep.depth_h>0.5 && pep.depth_h<1.5 && minimize<min3dG && minimize>min2dG){
            min3dG=minimize;
            third_fittest=pep;
        }

        if(pep.depth_h>0.5 && pep.depth_h<1.5 && minimize<min4dG && minimize>min3dG){
            min4dG=minimize;
            fourth_fittest=pep;
        }*/

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
    int cycle=0;
    while(cycle < 12) {  // Mutate sequences until 100 loops happen without finding a better fit

        Peptide newpep;
        Peptide temppep;

        cout<<fittest_pep.name<<" "<<second_fittest.name<<" "<<third_fittest.name<<" "<<fourth_fittest.name<<endl;

        crossover(fittest_pep,fourth_fittest, temppep.population);
        crossover(second_fittest,third_fittest, temppep.population);
        temppep.population.push_back(fittest_pep.name);
        temppep.population.push_back(second_fittest.name);
        temppep.population.push_back(third_fittest.name);
        temppep.population.push_back(fourth_fittest.name);
        int start2=temppep.population.size();
        maxdG=-9999999;
        max2dG=-9999999;
        max3dG=-9999999;
        max4dG=-9999999;

        for (int k=0;k<start2;++k){

            if (find_overlapping_SB(temppep.population[k])){
                Peptide mutations2;
                Peptide fittest_mut;
                mutation_3SB_avoid_overlap(mutations2.population, temppep.population[k], resids);
                /*int num=mutations.population.size();
                Peptide mutations2;
                mutations2.population.push_back(fittest_pep.name);
                mutations2.population.push_back(second_fittest.name);
                mutations2.population.push_back(third_fittest.name);
                mutations2.population.push_back(fourth_fittest.name);
                for (int m=0;m<num;m++){
                        add_mut(mutations2.population, mutations.population[m], resids);
                }*/

                double mx=-999999;
                for (int i=0;i<mutations2.population.size();i++) {

                    mutations2.name=mutations2.population[i];
                    mutations2.seq.resize(20);
                    mutations2.initial_pos("init_pos20", resids, res_id);

                    if(salt_bridges(mutations2) || aromatics(mutations2)){

                        vector<Peptide>SB_all;

                        population_SB(  mutations2, SB_all);
                        calc.get_energy_SB(SB_all,resids,res_id, mutations2);

                    }else{

                        calc.get_energy_total(mutations2, resids, res_id);

                    }

                        double maximize=mutations2.energy_h-mutations2.energy_b;
                        double score=maximize;

                        if (score>mx){
                            mx=score;
                            fittest_mut=mutations2;
                        }

                }
                newpep.population.push_back(fittest_mut.name);
            }else {
                Peptide mutations2;
                Peptide fittest_mut;
                mutation_NR(mutations2.population, temppep.population[k], resids);
                /*int num=mutations.population.size();
                Peptide mutations2;
                mutations2.population.push_back(fittest_pep.name);
                mutations2.population.push_back(second_fittest.name);
                mutations2.population.push_back(third_fittest.name);
                mutations2.population.push_back(fourth_fittest.name);
                for (int m=0;m<num;m++){
                    add_mut(mutations2.population, mutations.population[m], resids);
                }*/

                double mx=-999999;
                for (int i=0;i<mutations2.population.size();i++) {

                    mutations2.name=mutations2.population[i];
                    mutations2.seq.resize(20);
                    mutations2.initial_pos("init_pos20", resids, res_id);

                    if(salt_bridges(mutations2) || aromatics(mutations2)){

                        vector<Peptide>SB_all;

                        population_SB(  mutations2, SB_all);
                        calc.get_energy_SB(SB_all,resids,res_id, mutations2);

                    }else{

                        calc.get_energy_total(mutations2, resids, res_id);

                    }

                        double maximize=mutations2.energy_h-mutations2.energy_b;
                        double score=maximize;

                        if (score>mx){
                            mx=score;
                            fittest_mut=mutations2;
                        }

                }
                newpep.population.push_back(fittest_mut.name);

            }
        }

        Peptide one, two, three, four;
        int newstart=newpep.population.size();
        for (int i=0;i<newstart;i++) {
            //newpep.seq.clear();

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
                //double minimize=newpep.energy_h;
                double maximize=newpep.energy_h-newpep.energy_b;
                double max_z=newpep.depth_h-newpep.depth_b;
                double score=maximize;
                //double score=1000*max_z+0.1*maximize;
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
                    four=three;
                    three=two;
                    two=one;
                    one=newpep;
                }

                if(score>max2dG && score<maxdG){
                    max2dG=score;
                    four=three;
                    three=two;
                    two=newpep;
                }

                if(score>max3dG && score<max2dG){
                    max3dG=score;
                    four=three;
                    three=newpep;
                }

                if(score>max4dG && score<max3dG){
                    max4dG=score;
                    four=newpep;
                }
        /*if (newpep.depth_h>0.5 && newpep.depth_h<1.5 && minimize<mindG){
            mindG=minimize;
            fittest_pep=newpep;
        }

        if(newpep.depth_h>0.5 && newpep.depth_h<1.5 && minimize<min2dG && minimize>mindG){
            min2dG=minimize;
            second_fittest=newpep;
        }

        if(newpep.depth_h>0.5 && newpep.depth_h<1.5 && minimize<min3dG && minimize>min2dG){
            min3dG=minimize;
            third_fittest=newpep;
        }

        if(newpep.depth_h>0.5 && newpep.depth_h<1.5 && minimize<min4dG && minimize>min3dG){
            min4dG=minimize;
            fourth_fittest=newpep;
        }*/

            //cout<<newpep.name<<" "<<newpep.energy_h<<" "<<newpep.energy_b<<" " <<newpep.energy_h-newpep.energy_b<<endl;

        }
        fittest_pep=one;
        second_fittest=two;
        third_fittest=three;
        fourth_fittest=four;

        string peptide0=fittest_pep.name;

        // Sequence comparison with the last round
        if(peptide0==peptide){
            rep++;
        }else{
            rep=0;
        }
        cout<<fittest_pep.name<<" "<<rep<<" "<<fittest_pep.energy_h<<" "<<fittest_pep.energy_b<<" "<<fittest_pep.energy_h-fittest_pep.energy_b<<" "<<fittest_pep.depth_h<<" "<<fittest_pep.depth_b<<endl;

        if(rep==50){
            file1<<fittest_pep.name<<" "<<fittest_pep.energy_h<<" "<<fittest_pep.energy_b<<" "<<fittest_pep.energy_h-fittest_pep.energy_b<<" "<<fittest_pep.depth_h<<" "<<fittest_pep.depth_b<<endl;
            file1<<second_fittest.name<<" "<<second_fittest.energy_h<<" "<<second_fittest.energy_b<<" "<<second_fittest.energy_h-second_fittest.energy_b<<" "<<second_fittest.depth_h<<" "<<second_fittest.depth_b<<endl;
            file1<<third_fittest.name<<" "<<third_fittest.energy_h<<" "<<third_fittest.energy_b<<" "<<third_fittest.energy_h-third_fittest.energy_b<<" "<<third_fittest.depth_h<<" "<<third_fittest.depth_b<<endl;
            file1<<fourth_fittest.name<<" "<<fourth_fittest.energy_h<<" "<<fourth_fittest.energy_b<<" "<<fourth_fittest.energy_h-fourth_fittest.energy_b<<" "<<fourth_fittest.depth_h<<" "<<fourth_fittest.depth_b<<endl;

            Peptide mut_four;
            add_mut(mut_four.population, fittest_pep.name, resids);
            add_mut(mut_four.population, second_fittest.name, resids);
            add_mut(mut_four.population, third_fittest.name, resids);
            add_mut(mut_four.population, fourth_fittest.name, resids);

            fittest_pep.name=mut_four.population[0];
            second_fittest.name=mut_four.population[1];
            third_fittest.name=mut_four.population[2];
            fourth_fittest.name=mut_four.population[3];


            /*fittest_pep=second_fittest;
            second_fittest=third_fittest;
            third_fittest=fourth_fittest;
            fourth_fittest.name=gen_random();*/
            cycle=cycle+1;
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
