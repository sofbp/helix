#ifndef RESIDUE_H
#define RESIDUE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "atom.h"
#include <cassert>


using namespace std;

/*class Efile{
public:

    double z,E;

    Efile() {}
    Efile(double z, double E): z(z), E(E) {}

    // convert to array E[number_of_data_points].
    // Assumption of data points equidistant distribution -> fce( index ) -> z
};*/

class Residue
{
public:
    Residue() {
        Eh.resize(200);
        Eb.resize(200);
    }

    Atom pos;
    Atom CA;
    Atom com;
    //vector<Efile>Eh, Eb;

    double offset_Eh, offset_Eb;
    double factor_Eh, factor_Eb;
    double inv_factor_Eh, inv_factor_Eb;
    double upper_bound;
    double lower_bound;
    vector<double> Eh;
    vector<double> Eb;
    double cte=1.0/(0.008314*310);
    double e=2.718281828459045;

    void init(vector< vector<double>> Ehh, int type)
    {

        if(type==1){
            lower_bound=Ehh[199][0];
            upper_bound=Ehh[0][0];
            factor_Eh=(upper_bound-lower_bound)/199;
            inv_factor_Eh=1.0/factor_Eh;
            for (int i=0; i<200; ++i) {
                this->Eh[i] = Ehh[i][1];
            }
        }else{
            lower_bound=Ehh[199][0];
            upper_bound=Ehh[0][0];
            factor_Eb=(upper_bound-lower_bound)/199;
            inv_factor_Eb=1.0/factor_Eb;
            for (int i=0; i<200; ++i) {
                this->Eb[i] = Ehh[i][1];
            }
        }
    }


    void E_charged_init( Residue b, vector<vector<double>>& newEh, int type)
    {
        if(type==1){
            double newE, Gch, Gn;
            for (int i=0;i<200;i++) {
                double z=(199-i)*factor_Eh;
                if(z>upper_bound){
                    newE=0;
                }else{
                    int convert0 = int ((upper_bound-z) * inv_factor_Eh);
                    int convert1 = int ((b.upper_bound-z) * b.inv_factor_Eh);


                    //double proportional_remainder0 = fmod(z, factor_Eh) * inv_factor_Eh;
                    //double proportional_remainder0 = (z-(int(z*inv_factor_Eb) * factor_Eb))*inv_factor_Eb;
                    double proportional_remainder0 = (z * inv_factor_Eh) - int(z *inv_factor_Eh);
                    //double proportional_remainder1 = fmod(z, b.factor_Eh) * b.inv_factor_Eh;
                    //double proportional_remainder1 = (z-(int(z*b.inv_factor_Eb) * b.factor_Eb))*b.inv_factor_Eb;
                    double proportional_remainder1 = (z * inv_factor_Eh) - int(z *inv_factor_Eh);

                    Gch=(1.0-proportional_remainder1)*b.Eh[convert1] + proportional_remainder1*b.Eh[convert1+1];
                    Gn=(1.0-proportional_remainder0)*Eh[convert0] + proportional_remainder0*Eh[convert0+1];

                    double total=pow(e, -Gn*cte)+pow(e, -Gch*cte);
                    double prop=pow(e, -Gn*cte)/total;

                    newE = (1.0-prop)*Gch + prop*Gn;
                }
                newEh.resize(200, vector<double> (2));
                newEh[i][0]=z;
                newEh[i][1]=newE;
                //cout<<z<<" "<<Gch<<" "<<Gn<<" "<<newE<<endl;
            }
        }else {
            double newE, Gch, Gn;
            for (int i=0;i<200;i++) {
                double z=(199-i)*factor_Eb;
                if(z>upper_bound){
                    newE=0;
                }else{
                    int convert0 = int ((upper_bound-z) * inv_factor_Eb);
                    int convert1 = int ((b.upper_bound-z) * b.inv_factor_Eb);


                    //double proportional_remainder0 = fmod(z, factor_Eb) * inv_factor_Eb;
                    double proportional_remainder0 = (z * inv_factor_Eb) - int(z *inv_factor_Eb);
                    //double proportional_remainder1 = fmod(z, b.factor_Eb) * b.inv_factor_Eb;
                    //double proportional_remainder1 = (z-(int(z*b.inv_factor_Eb) * b.factor_Eb))*b.inv_factor_Eb;
                    double proportional_remainder1 = (z * inv_factor_Eb) - int(z *inv_factor_Eb);

                    Gch=(1.0-proportional_remainder1)*b.Eb[convert1] + proportional_remainder1*b.Eb[convert1+1];
                    Gn=(1.0-proportional_remainder0)*Eb[convert0] + proportional_remainder0*Eb[convert0+1];

                    double total=pow(e, -Gn*cte)+pow(e, -Gch*cte);
                    double prop=pow(e, -Gn*cte)/total;

                    newE = (1.0-prop)*Gch + prop*Gn;
                }
                newEh.resize(200, vector<double> (2));
                newEh[i][0]=z;
                newEh[i][1]=newE;
                //cout<<z<<" "<<Gch<<" "<<Gn<<" "<<newE<<endl;
            }

        }
    }

    inline double get_E_human(double z)
    {

        if(z >= upper_bound) {
            return 0.0;
        }

        if(z < lower_bound) {
            z=z*-1;
        }


        int convert = int ((upper_bound-z) * inv_factor_Eh);
        //cout<<upper_bound<<" "<<lower_bound<< " "<<convert<<" "<<z<< " "<< int ((upper_bound-0) * inv_factor_Eh)<<endl;
        //assert(convert <= 198);
        //double proportional_remainder = (z-(int(z*inv_factor_Eh) * factor_Eh))*inv_factor_Eh;
        double proportional_remainder = (z * inv_factor_Eh) - int(z *inv_factor_Eh);
        //double proportional_remainder = fmod(z, factor_Eh) * inv_factor_Eh;
        //cout<<" "<<new_rem<<" "<<proportional_remainder<<endl;
        //exit(0);
        if(convert >= 199){
            //cout<<"here"<<endl;
            return (Eh[convert]);
        }else{
            //assert(convert <= 198);
            return  (1.0-proportional_remainder)*Eh[convert] + proportional_remainder*Eh[convert+1];

        }
    }

    inline double get_E_bacteria(double z)
    {
        if(z >= upper_bound) {
            return 0.0;
        }

        if(z < lower_bound) {
            z=z*-1;
        }

        int convert = int ((upper_bound-z) * inv_factor_Eb);
        //assert(convert <= 198);
        //double proportional_remainder = (z-(int(z*inv_factor_Eb) * factor_Eb))*inv_factor_Eb;
        double proportional_remainder = (z * inv_factor_Eb) - int(z *inv_factor_Eb);
        //double proportional_remainder = fmod(z, factor_Eb) * inv_factor_Eb;
        if(convert >= 199){
            //cout<<"here"<<endl;
            return (Eb[convert]);
        }else{
            //assert(convert <= 198);
            return  (1.0-proportional_remainder)*Eb[convert] + proportional_remainder*Eb[convert+1];

        }
    }

    inline double get_E_human(double z, vector< vector<double>>& Eh)
    {

        lower_bound=Eh[199][0];
        upper_bound=Eh[0][0];
        factor_Eh=(upper_bound-lower_bound)/199;

        if(z >= upper_bound) {
            return 0.0;
        }

        if(z < lower_bound) {
            z=z*-1;
        }

        int convert = int ((upper_bound-z)/factor_Eh);
        double proportional_remainder = fmod(z, factor_Eh) / factor_Eh;
        //cout<<convert<<" "<<z<<" "<<Eh[convert][0]<<" "<<Eh[convert+1][0]<<endl;
        return  (1.0-proportional_remainder)*Eh[convert][1] + proportional_remainder*Eh[convert+1][1];
    }
};


#endif // RESIDUE_H
