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

using namespace std;

const double radToDeg = 57.2958;
const double degToRad = 1.0/57.2958;

class Efile{
public:

    double z,E;

    Efile() {}
    Efile(double z, double E): z(z), E(E) {}
};



void load_sequnce(string input, vector<string>& population)
{
    std::fstream fs( input, std::fstream::in );

    string c1;
    while( !fs.eof() ) // Lines in input, eof = end of file
    {
        fs >> c1 ;
        population.push_back(c1);
    }
    population.pop_back();
}

class calculate_ddG{
public:
    calculate_ddG() {}



    void initial_pos(string input, vector<Atom>& seq) // sequence manipulation
    {
        std::fstream init( input, std::fstream::in );
        double xpos,ypos,zpos;

        while( !init.eof() ) // Lines in input, For (steps), eof = end of file
        {
            Atom add;
            init >> xpos >> ypos >> zpos;
            add.x=xpos;
            add.y=ypos;
            add.z=zpos;
            seq.push_back(add);
        }
        seq.pop_back();

    }



    void get_aa(vector<string> population, int i, vector<string>& aa){
        std::stringstream pop ( population[i] );
        char c;
        string c1;
        while( pop.get(c) )
            {
                c1=c;
                    aa.push_back(c1);
            }
    }

    void energy_file(vector<Efile>& NW, string aa, string tp){

        string nw=aa+tp;

        std::fstream fs( nw, std::fstream::in );
        double xx,yy;
        Efile add;

        if(!fs.good()){
            cout<<nw<<" Unknown"<<endl;
            exit(1);
        }

        while( !fs.eof() ) // Lines in input, eof = end of file
        {
            fs >> xx >> yy;
            add.z=xx;
            add.E=yy;
            NW.push_back(add);
        }
        NW.pop_back();
    }

    double get_energy(int i, double x_pos, vector<Efile> En){

        double post, prev;
        double G=0.00;
        int start=i*200;
        int end=start+200;


        if(x_pos<0){
            x_pos=x_pos*-1;
        }

        for (int n=start;n<end;n++){
            if(x_pos < En[n].z && x_pos > En[n+1].z){
                post = En[n].E;
                prev = En[n+1].E;
                G = (prev+post)/2;
            }
        }
        return G;
    }

    void get_energy_total(int rot_deg, int tilt_deg, int depth, vector<Atom> seq, vector<double>& energy, vector<Efile> E){
        Atom current;
        Atom axis_y=Atom(0,1,0);
        Atom axis_x=Atom(1,0,0);
        axis_x.normalise();
        axis_y.normalise();
        double total_en,en=0;

        for (int i=0;i<rot_deg;i=i+20) {
            double deg=i*degToRad;
            for (int j=0;j<tilt_deg;j=j+20) {
                double deg1=j*degToRad;
                //printf ("%d%03d ", j, i);
                //cout<<i<<j<<" ";
                for (int k=0;k<depth;++k) {
                    double disp=k*0.1;
                     total_en=0;
                    for (int a=0;a<seq.size();++a) { //move all beads

                        current=seq[a];

                        //Rotate, tilt, translate -->traslations must be the last because rotations are done relative to 0,0,0.

                        current.rotate(axis_y, deg);
                        current.rotate(axis_x, deg1);
                        current.z=current.z+disp;

                        //calculate deltaG for new positions
                        en=get_energy(a, current.z, E);
                        total_en+=en;
                    }
                    energy.push_back(total_en);
                }
            }
        }
    }


    double deltaG(vector<double>energy, int pos, int pos2, int wat){
        double e1,e2,dG;
        int idx1,idx2;
        double dGbig=-1;

        for (int p=pos;p<pos2;p++) {

            for (int i=0;i<energy.size();++i) {
                idx1=i*wat+wat;
                idx2=i*p+p;

                if(idx1>energy.size() || idx2 > energy.size()){
                    break;
                }
                e1=energy[idx1];
                e2=energy[idx2];
                dG= e2-e1;

                if(dG>0){
                    //cout<<i%rot_deg*20%rot_deg<<" "<<i/10*20%rot_deg<<" "<<dG<<endl;  //print angles for each energy
                    if(dG>dGbig){
                        dGbig=dG;
                    }
                }
            }
        }
        return dGbig;
    }

    double deltaGlow(vector<double>energy, int pos, int pos2, int wat){
        double e1,e2,dG;
        int idx1,idx2;
        double dGsmall=0;

        for (int p=pos;p<pos2;p++) {

            for (int i=0;i<energy.size();++i) {
                idx1=i*wat+wat;
                idx2=i*p+p;

                if(idx1>energy.size() || idx2 > energy.size()){
                    break;
                }
                e1=energy[idx1];
                e2=energy[idx2];
                dG= e2-e1;

                if(dG<0){
                    //cout<<i%rot_deg*20%rot_deg<<" "<<i/10*20%rot_deg<<" "<<dG<<endl;  //print angles for each energy
                    if(dG<dGsmall){
                        dGsmall=dG;
                    }
                }
            }
        }
        return dGsmall;
    }

    void crossover(vector<string>& pep1, vector<string>& pep2){
        //srand (time(NULL));
        int cross_point=rand() %20;
        for (int i = 0; i < cross_point; i++) {
            string temp = pep1[i];
            pep1[i] = pep2[i];
            pep2[i] = temp;

        }
    }

    void analyze_G(vector<string>population, vector<string>& pep){

        pep.clear();
        //double ddG;
        double small=0;
        double small2=0;
        vector<string> pep1, pep2;

        //Get alpha carbon positions from file
        vector<Atom>seq;
        initial_pos("init_pos20", seq);

        for (int i=0;i<population.size();++i) {
            vector<string> aa;
            get_aa(population, i, aa);

            if(aa.size() != seq.size()){
                cout<< "Input files have different length."<<endl;
                exit(0);
            }

            Atom CM=com(seq);
            vector<Efile>Eh, Eb;

            //move peptide COM to (0,0,0) --> starting point
            for (int j=0;j<seq.size();j++) {
                seq[j]=seq[j]-CM;
                energy_file(Eh,aa[j], "h");
                //energy_file(Eb,aa[j], "b");
            }

            int rot_deg=200;
            int tilt_deg=200;
            int depth=30;

            vector<double>energy_h, energy_b;
            //auto t1 = std::chrono::high_resolution_clock::now();

            get_energy_total(rot_deg, tilt_deg, depth, seq, energy_h, Eh);
            //get_energy_total(rot_deg, tilt_deg, depth, seq, energy_b, Eb);

            int wat=30;
            int pos=15;
            int pos2=20;
            double dGh, dGb;

            //Energy of adsorption. Take min energy for bacterial membrane, and max energy for human membrane.
            dGh=deltaG(energy_h, pos, pos2, wat);
            //dGb=deltaGlow(energy_b, pos, pos2, wat);
            //ddG=dGh-dGb;


            if(dGh>small){
                small=dGh;
                pep1=aa;
            }
            if(dGh>small2 && dGh < small){
                small2=dGh;
                pep2=aa;
            }

            //cout<<dGh<<endl;
            //print_map(energy_h, depth,rot_deg);

            //auto t2 = std::chrono::high_resolution_clock::now();

            //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

            /*std::cout << duration<<endl;
        cout<<"done"<<endl;*/
        }
        //crossover
        /*crossover(pep1,pep2);*/

        //mutation

        pep=pep1;

        //cout<<small<<" "<<small2<<endl;
    }



};

string print_peptide(vector<string>pep){
    string peptide;
    for (int n=0;n<pep.size();n++) {
        peptide=peptide+pep[n];
    }
    return peptide;
}

void mutation(vector<string>& newpopulation,  vector<string> pep){

    vector<string> resids = {"A", "E", "K","W", "R", "C", "V", "L", "I", "M", "Q", "F", "N"};
    //srand (time(NULL));
    int mut_point=rand() %20;

    for (int i=0;i<resids.size();i++) {
        pep[mut_point]=resids[i];
        string newpeptide=print_peptide(pep);
        newpopulation.push_back(newpeptide);
    }
}

int main()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    calculate_ddG ddG;
    vector<string> population,pep;

    //Get amino acid sequence from input file

    load_sequnce("input", population);

    //Initial population analysis --> select best fit
    ddG.analyze_G(population, pep);
    string peptide0=print_peptide(pep);

    //Mutation of best fit --> new set of peptides --> analysis and best fit selection
    int rep=0;
    while(rep<30) {

        vector<string> newpopulation;

        mutation(newpopulation,pep);

        ddG.analyze_G(newpopulation, pep);
        string peptide=print_peptide(pep);

        if(peptide0==peptide){
            rep++;
        }else{
            rep=0;
        }
        cout<<peptide<<" "<<rep<<endl;
        peptide0=peptide;
        peptide.clear();

    }
    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    std::cout << duration<<endl;

    return 0;
}
