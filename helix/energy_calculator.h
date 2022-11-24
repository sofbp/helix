#ifndef ENERGY_CALCULATOR_H
#define ENERGY_CALCULATOR_H

#include <vector>
#include "atom.h"
#include "peptide.h"
#include "residue.h"
#include "chrono"
#include <sys/time.h>
#include "map"

/*
 * std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
 * std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
 * std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
*/

using namespace std::chrono;

class Ener
{
public:
    double energy;
    double tilt;
};

class Energy_calculator
{
private:
    const double radToDeg = 57.2958;
    const double degToRad = 1.0/57.2958;

public:
    Energy_calculator() {}

    int rot_deg=360;
    int tilt_deg=360;
    int depth=40;
    vector< vector< vector< double > > > E_H;
    vector< vector< vector< double > > > E_B;
    vector< vector< vector< double > > > Charged;
    vector< vector< vector< double > > > Charged_B;
    double cte=1.0/(0.008314*310);
    double e=2.718281828459045;


    void load_EfilesH(char input, int n)
    {
        string in;
        in.push_back(input);
        string human=in+"h";
        std::fstream fs( human, std::fstream::in );
        double zz, ee;

        int i=0;
        while( !fs.eof() && i<200 ) // Lines in input, eof = end of file
        {
            fs >> zz >> ee;
            E_H[n][i][0]=zz;
            E_H[n][i][1]=ee;

            ++i;
        }
    }

    void load_EfilesB(char input, int n)
    {
        string in;
        in.push_back(input);
        string bact=in+"b";
        std::fstream fs( bact, std::fstream::in );
        double zz, ee;

        int i=0;
        while( !fs.eof() && i<200 ) // Lines in input, eof = end of file
        {
            fs >> zz >> ee;
            E_B[n][i][0]=zz;
            E_B[n][i][1]=ee;

            ++i;
        }
    }

    void load_E_charged(string input, int n)
    {

        std::fstream fs( input, std::fstream::in );
        double zz, ee;
        int i=0;
        while( !fs.eof() && i<200 ) // Lines in input, eof = end of file
        {

            fs >> zz >> ee;
            Charged[n][i][0]=zz;
            Charged[n][i][1]=ee;
            ++i;
        }
    }

    void init_grid(int a, int b, int c)
    {
        E_H.resize(a, vector<vector<double>>(b, vector<double>(c)));
        for (int x=0;x<a;++x) {
            for (int y=0;y<b;++y) {
                for (int z=0;z<c;++z) {
                    E_H[x][y][z]=0;

                }
            }
        }
        E_B.resize(a, vector<vector<double>>(b, vector<double>(c)));
        for (int x=0;x<a;++x) {
            for (int y=0;y<b;++y) {
                for (int z=0;z<c;++z) {
                    E_B[x][y][z]=0;

                }
            }
        }

    }

    void get_energy_total(Peptide& current, Peptide aa, map<char, int> res_id){
        int interval=5;
        Atom axis_y=Atom(0,1,0);
        Atom axis_x=Atom(1,0,0);
        vector<vector<double>>energies_h;
        energies_h.resize(depth, vector<double> (rot_deg*tilt_deg/interval));
        vector<vector<double>>energies_b;
        energies_b.resize(depth, vector<double> (rot_deg*tilt_deg/interval));
        vector<double>G_all;
        vector<double>G_all_b;
        //cout<<energies_h[0].size()<<endl;
        //axis_x.normalise();
        //axis_y.normalise();
        double min=999999;
        double minB=999999;
        double total_B_en_h=0;
        double total_B_en_b=0;
        double prob_sum=0;
        double prob_sumB=0;
        double shift=0;
        double shiftB=0;

        double total_en_h;
        double total_en_b;
        for (int k=0;k<depth;++k) {
            double disp=k*0.1;
            for (int i=0;i<rot_deg;i=i+interval) {
                double deg_rot=i*degToRad;

                for (int j=0;j<tilt_deg;j=j+interval) {
                    double deg_tilt=j*degToRad;
                    //printf ("%d%03d ", j, i);

                    total_en_h=0;
                    total_en_b=0;

                    //
                    // Move, rotate and get E under peptide class
                    //
                    for (int a=0;a<current.size();++a) {

                        Atom positions;
                        positions.x=current.seq[a].pos.x;
                        positions.y=current.seq[a].pos.y;
                        positions.z=current.seq[a].pos.z;

                        //Rotate, tilt, translate -->traslations must be the last because rotations are done relative to 0,0,0.

                        positions.rotate(axis_y, deg_rot);
                        positions.rotate(axis_x, deg_tilt);
                        positions.z=positions.z+disp;

                        //calculate deltaG for new positions

                        int inde=res_id[current.name[a]];
                        //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


                        double en_h=aa.seq[inde].get_E_human(positions.z);

                        double en_b=aa.seq[inde].get_E_bacteria(positions.z);
                        //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;

                        total_en_h+=en_h;
                        total_en_b+=en_b;
                    }
                    //
                    // Boltzmann average --> average out degrees of freedom from tilt and rotation
                    // Calculate P(x,y) = e^(-dG/kT)
                    //
                    double B_en_h=pow(e,-total_en_h*cte);
                    energies_h[k][(i*j+j)/interval]=B_en_h;
                    //energies_h[k].push_back(B_en_h);
                    //cout<<k<<" "<<(i*j+j)/interval<<endl;
                    total_B_en_h+=B_en_h;

                    double B_en_b=pow(e,-total_en_b*cte);
                    //energies_b[k].push_back(B_en_b);
                    energies_b[k][(i*j+j)/interval]=B_en_b;
                    total_B_en_b+=B_en_b;
                }
                //exit(1);
            }
        }
        //cout<<energies_h[0].size()<<endl;

        current.totalP_h=total_B_en_h;
        current.totalP_b=total_B_en_b;

        //
        // P(x) = sum (P(x,y)/total_sum) --> For each z value
        // G(x) = -kT*ln(P(x))
        //
        double G, Gb;
        for (int q=0;q<depth;q++) {

            for (int p=0;p<energies_h[q].size();p++) {

                double prob=energies_h[q][p]/total_B_en_h;
                double probB=energies_b[q][p]/total_B_en_b;

                prob_sumB+=probB;

                prob_sum+=prob;

            }
            G=-cte*log(prob_sum);
            G_all.push_back(G);
            Gb=-cte*log(prob_sumB);
            G_all_b.push_back(Gb);

            if(q==(depth-1)){
                shift=G;
                shiftB=Gb;
            }
            prob_sum=0;
            prob_sumB=0;

        }

        //
        // Shift G profile so that it is 0 at z=4 (outside the membrane)
        // Store energy minima and depth
        //

        vector<double> newG, newGb;
        double z, zb;

        for (int r=0;r<depth;r++) {
            newG.push_back(G_all[r]-shift);
            newGb.push_back(G_all_b[r]-shiftB);
            z=r*0.1;
            zb=r*0.1;
            if(newG[r]<min){
                min=newG[r];

                current.energy_h=min;
                current.depth_h=z;
            }

            if(newGb[r]<minB){
                minB=newGb[r];

                current.energy_b=minB;
                current.depth_b=zb;

            }
        }
    }

    void print_Gplot(Peptide& current, Peptide aa, map<char, int> res_id, string name){
        Atom axis_y=Atom(0,1,0);
        Atom axis_x=Atom(1,0,0);
        vector<vector<double>>energies_h;
        energies_h.resize(depth, vector<double> (rot_deg*tilt_deg/5));
        vector<vector<double>>energies_b;
        energies_b.resize(depth, vector<double> (rot_deg*tilt_deg/5));
        vector<double>G_all;
        vector<double>G_allB;
        //axis_x.normalise();
        //axis_y.normalise();
        double min=999999;
        double total_B_en=0;
        double total_B_en_b=0;
        double prob_sum=0;
        double prob_sumB=0;
        double shift=0;
        double shiftB=0;

        string fileH=name+"_GH";
        string fileB=name+"_GB";

        double total_en_h=0;
        double total_en_b=0;
        for (int k=0;k<depth;++k) {
            double disp=k*0.1;
            for (int i=0;i<rot_deg;i=i+5) {
                double deg_rot=i*degToRad;

                for (int j=0;j<tilt_deg;j=j+5) {
                    double deg_tilt=j*degToRad;
                    //printf ("%d%03d ", j, i);

                    total_en_h=0;
                    total_en_b=0;

                    //
                    // Move rotates and get E under peptide class
                    //
                    for (int a=0;a<current.size();++a) { //move all beads

                        Atom positions;
                        positions.x=current.seq[a].pos.x;
                        positions.y=current.seq[a].pos.y;
                        positions.z=current.seq[a].pos.z;

                        //Rotate, tilt, translate -->traslations must be the last because rotations are done relative to 0,0,0.

                        positions.rotate(axis_y, deg_rot);
                        positions.rotate(axis_x, deg_tilt);
                        positions.z=positions.z+disp;

                        //calculate deltaG for new positions
                        //double en=current.seq[a].get_E_human(positions.z, E[res_id[current.seq[a].type]]);

                        int inde=res_id[current.name[a]];
                        //cout<<i<<" "<<j<<" "<<k<<" "<<current.name[a]<<endl;

                        double en_h=aa.seq[inde].get_E_human(positions.z);
                        double en_b=aa.seq[inde].get_E_bacteria(positions.z);

                        total_en_h+=en_h;
                        total_en_b+=en_b;
                    }
                    double B_en=pow(e,-total_en_h*cte);
                    energies_h[k].push_back(B_en);

                    double B_en_b=pow(e,-total_en_b*cte);
                    energies_b[k].push_back(B_en_b);

                    total_B_en+=B_en;
                    total_B_en_b+=B_en_b;

                }
            }
        }


        double G, Gb;
        for (int q=0;q<depth;q++) {

            for (int p=0;p<energies_h[q].size();p++) {

                double prob=energies_h[q][p]/total_B_en;
                double probB=energies_b[q][p]/total_B_en_b;

                prob_sumB+=probB;

                prob_sum+=prob;

            }
            G=-cte*log(prob_sum);
            G_all.push_back(G);
            Gb=-cte*log(prob_sumB);
            G_allB.push_back(Gb);

            if(q==(depth-1)){
                shift=G;
                shiftB=Gb;
            }
            prob_sum=0;
            prob_sumB=0;

        }

        double z;
        ofstream myfile2;
        myfile2.open (fileH);
        ofstream myfile3;
        myfile3.open (fileB);

        for (int r=0;r<depth;r++) {
            z=r*0.1;
            myfile2<<z<<" "<<(G_all[r]-shift)<<endl;
            myfile3<<z<<" "<<(G_allB[r]-shiftB)<<endl;

        }
        myfile2.close();
        myfile3.close();
    }

    void print_Emap(Peptide& current, Peptide aa, map<char, int> res_id, string name){
        Atom axis_y=Atom(0,1,0);
        Atom axis_x=Atom(1,0,0);
        vector<vector<double>>energies_h;
        energies_h.resize(depth, vector<double> (rot_deg*tilt_deg/5));
        vector<vector<double>>energies_b;
        energies_b.resize(depth, vector<double> (rot_deg*tilt_deg/5));
        vector<double>G_all;
        //axis_x.normalise();
        //axis_y.normalise();
        double min=999999;
        double total_B_en=0;
        double total_B_en_b=0;
        ofstream myfile, myfile2,myfile3, myfile4;
        string fileH = name + "_H";
        string fileB = name + "_B";
        myfile.open (fileH);
        //myfile2.open ("P_mapH");
        myfile3.open (fileB);
        //myfile4.open ("P_mapB");

        myfile<<depth<<" ";
        //myfile2<<depth<<" ";
        myfile3<<depth<<" ";
        //myfile4<<depth<<" ";

        for (int n=0;n<depth;++n) {
            myfile<<n<<" ";
            //myfile2<<n<<" ";
            myfile3<<n<<" ";
            //myfile4<<n<<" ";
        }
        myfile<<endl;
        //myfile2<<endl;
        myfile3<<endl;
        //myfile4<<endl;
        double total_en_h=0;
        double total_en_b=0;
            for (int i=0;i<rot_deg;i=i+5) {
                double deg_rot=i*degToRad;

                for (int j=0;j<tilt_deg;j=j+5) {
                    double deg_tilt=j*degToRad;
                    //printf ("%d%03d ", j, i);

                    myfile<<j<<i<<" ";
                    myfile2<<j<<i<<" ";
                    myfile3<<j<<i<<" ";
                    myfile4<<j<<i<<" ";
                    for (int k=0;k<depth;++k) {
                        double disp=k*0.1;

                    total_en_h=0;
                    total_en_b=0;

                    //
                    // Move rotates and get E under peptide class
                    //
                    for (int a=0;a<current.size();++a) { //move all beads

                        Atom positions;
                        positions.x=current.seq[a].pos.x;
                        positions.y=current.seq[a].pos.y;
                        positions.z=current.seq[a].pos.z;

                        //Rotate, tilt, translate -->traslations must be the last because rotations are done relative to 0,0,0.

                        positions.rotate(axis_y, deg_rot);
                        positions.rotate(axis_x, deg_tilt);
                        positions.z=positions.z+disp;

                        //calculate deltaG for new positions
                        //double en=current.seq[a].get_E_human(positions.z, E[res_id[current.seq[a].type]]);

                        int inde=res_id[current.name[a]];
                        //cout<<i<<" "<<j<<" "<<k<<" "<<current.name[a]<<endl;

                        double en_h=aa.seq[inde].get_E_human(positions.z);
                        double en_b=aa.seq[inde].get_E_bacteria(positions.z);

                        total_en_h+=en_h;
                        total_en_b+=en_b;
                    }
                    double B_en=pow(e,-total_en_h*cte);
                    energies_h[k].push_back(B_en);

                    double B_en_b=pow(e,-total_en_b*cte);
                    energies_b[k].push_back(B_en_b);

                    total_B_en+=B_en;
                    total_B_en_b+=B_en_b;
                    myfile<<total_en_h<<" ";
                   // myfile2<<B_en/current.totalP_h<<" ";
                    myfile3<<total_en_b<<" ";
                    //myfile4<<B_en_b/current.totalP_b<<" ";

                }
                    myfile<<endl;
                    //myfile2<<endl;
                    myfile3<<endl;
                   // myfile4<<endl;
            }
        }
            myfile.close();
            //myfile2.close();
            myfile3.close();
           // myfile4.close();

    }

};

#endif // ENERGY_CALCULATOR_H
