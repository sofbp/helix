#include <QCoreApplication>
#include <beads.h>
#include <string>
#include <iostream>
#include <cmath>

class Atom{
public:

    double x,y,z;

    string type;
    Atom() {}
    Atom(double x, double y, double z): x(x), y(y), z(z) {}
    Atom(double x, double y, double z, string type): x(x), y(y), z(z), type(type) {}

    bool operator==(const Atom& o) const {
        if(x != o.x) return false;
        if(y != o.y) return false;
        if(z != o.z) return false;

        return true;
    }

    bool operator!=(const Atom& o) const {
        return !(*this == o);
    }

    Atom operator+(const Atom& o) const {
        return Atom(x+o.x, y+o.y, z+o.z);
    }

    void operator+=(const Atom& o) {
        x += o.x;
        y += o.y;
        z += o.z;
    }

    Atom operator-(const Atom& o) const {
        return Atom(x-o.x, y-o.y, z-o.z);
    }

    Atom operator*(const double a) const {
        return Atom(x*a,y*a,z*a, type);
    }

    Atom operator/(const double a) const {
        return Atom(x/a,y/a,z/a, type);
    }

    void operator*=(double a) {
        x*=a;
        y*=a;
        z*=a;
    }

    inline void rotate(Atom& axis, double angle) {
        angle*=0.5;
        double cosAngle = cos(angle);
        double sinAngle = sin(angle);
        double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz;
        double qw = cosAngle, qx = (axis.x * sinAngle), qy = (axis.y * sinAngle), qz = (axis.z * sinAngle);

        /*    t1 = quat.w * quat.w; */
        t2 =  qw * qx;
        t3 =  qw * qy;
        t4 =  qw * qz;
        t5 = -qx * qx;
        t6 =  qx * qy;
        t7 =  qx * qz;
        t8 = -qy * qy;
        t9 =  qy * qz;
        t10 = -qz * qz;

        newx = 2.0 * ( (t8+t10) * x + (t6-t4)  * y + (t3+t7) * z ) + x;
        newy = 2.0 * ( (t4+t6)  * x + (t5+t10) * y + (t9-t2) * z ) + y;
        newz = 2.0 * ( (t7-t3)  * x + (t2+t9)  * y + (t5+t8) * z ) + z;

        x = newx;
        y = newy;
        z = newz;
    }

    int rot(int i){
        //rotate z, x axis
    }

    int tilt(int j){
        //rotate z,y axit
    }

};

class Input{
public:
    Input() {}

    string a1, a2, a3, a4;

    bool loadInput(string input)
    {
        cout<<"hola";
        std::fstream fs( input, std::fstream::in );
        string what;
        stringstream ss;

        while( !fs.eof() ) // Lines in input
        {
            ss.flush();
            ss.clear();

            ss >> what;
            if( what.compare("Seq:") == 0 )   { ss >> a1 >> a2 >> a3 >> a4 ; }

            what.clear();
        }
        fs.close();

        return true;
    }

};

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

int main(int argc, char *argv[])
{
    cout << argc << " " << argv[0] << " " << argv[1] << endl;;
    return 0;
    /*Input b;
    b.loadInput("input");*/
    vector<string> aa;
    //cout<<b.a1<<" " <<b.a2<<" "<<b.a3<<endl;
    aa.push_back("A");
    aa.push_back("L");
    //cout<<seq[0]<<endl;
    Atom first;
    int x=1, y=1, z=0;
    first.x=x;
    first.y=y;
    first.z=z;
    first.type=aa[0];
    vector<Atom>seq;
    seq.push_back(first); //seq sera el vector que contenga las caracteristicas de las beads(x y z type)

    for (int n=1; n<aa.size(); ++n){
        Atom prev = seq[n-1];
        Atom next;
        //Aqui deberia calcularse la posicion del siguiente aa, rotacion y traslacion, en funcion del anterior
        next.x=prev.x;
        next.y=prev.y+1.5;
        next.z=prev.z;
        seq.push_back(next);
    }

    //calc COM of seq, all movement in reference to that
    Atom CM=com(seq);

    int rot_deg=360;
    int tilt_deg=360;
    int depth=3;
    Atom axis_y=Atom(0,1,0);
    Atom axis_x=Atom(1,0,0);

    for (int i=0;i<rot_deg;++i) {
        for (int j=0;j<tilt_deg;++j) {
            for (int k=0;k<depth;++k) {

                for (int a=0;a<seq.size();++a) { //move all beads
                    Atom current=seq[a];
                    current.z=current.z+k;
                    current.rotate(axis_y, i); //pero quiero rotar sobre el com del peptido???????
                    current.rotate(axis_x,j);
                    //calculate deltaG for new positions
                    //store deltaG of peptide with its rot,tilt and tras values
                }
            }
        }
    }
    return 0;
}



