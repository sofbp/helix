#ifndef POPULATION_H
#define POPULATION_H

#include "peptide.h"
#include "energy_calculator.h"

using namespace std;

class Population
{
public:
    Population() {}

    double maxdG=-9999999;
    double max2dG=-9999999;
    Peptide fittest_pep, second_fittest;

    void population_analysis(Peptide& pep){
        Energy_calculator calc;

        for (int i=0;i<pep.population.size();i++) {
            pep.seq.clear();
            pep.get_aa(i);
            pep.initial_pos("init_pos20");
            calc.get_energy_total(pep);

            //cout<<pep.population[i]<<" "<< pep.energy<<" "<<pep.rot_angle<<" "<<pep.tilt_angle<<" "<<pep.depth*0.1<<endl;

            // selectivity criteria
            if(pep.energy>maxdG){
                maxdG=pep.energy;
                fittest_pep=pep;
            }
            if(pep.energy>max2dG && pep.energy < maxdG){
                max2dG=pep.energy;
                second_fittest=pep;
            }

            /*for (int i=0;i<calc.energy.size();i++) {
        cout<<pep.rot_angle[i]<<" "<<pep.tilt_angle[i]<<" "<<pep.depth[i]<<" "<<calc.energy[i]<<endl;
        //cout<<pep.rot_angle.size()<<endl;
    }*/
        }
    }
};

#endif // POPULATION_H
