#ifndef PEPTIDE_H
#define PEPTIDE_H

#include <vector>
#include "residue.h"

using namespace std;

class Peptide
{
public:
    Peptide() {}

    vector<Residue> pep;

    Atom com( vector<Atom> seq)
    {
        Atom cm;
        int count=0;
        for(int i=0; i<seq.size(); ++i)
        {
                    cm.x += seq[i].x;
                    cm.y += seq[i].y;
                    cm.z += seq[i].z;
                    ++count;
        }

        cm *= 1.0/count;
        return cm;
    }
};

#endif // PEPTIDE_H
