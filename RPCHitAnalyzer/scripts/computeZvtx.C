#include "Math/LorentzVector.h"
#include "Math/Vector4Dfwd.h"


class quadraticEquation
{
public:
  quadraticEquation(double vala, double valb, double valc)
    : a(vala), b(valb), c(valc)
  {
    Delta=b*b-4*a*c;
    f1=-b/(2*a);
    f2=0;
    if (hasSolution()) f2=sqrt(Delta)/(2*a);
  }
  
  bool hasSolution() {return Delta>=0;}
  bool oneSolution() {return Delta==0;}
  int Nsolution() { if (Delta >0) return 2; else if (Delta == 0) return 1; else return 0;}
  double solmoins() {return f1-f2;}
  double solplus() {return f1+f2;}
private:
  double a,b,c;
  double Delta;
  double f1,f2;
};

class ZTcalc
{
public:
  ZTcalc(ROOT::Math::XYZTVector a,ROOT::Math::XYZTVector b)
    : hita(a), hitb(b) {init();}

  
  ZTcalc(double xa,double ya, double za, double ta,
	 double xb,double yb, double zb, double tb,
	 double cspeed=29.9792458 /* cm/ns */)
    : hita(xa,ya,za,ta*cspeed), 
      hitb(xb,yb,zb,tb*cspeed)
  {init();}
  
  void compute();
  int Nsolution() {return Nsol;}
  double Zsol(int i) {return Z[i];}
  double Tsol(int i) {return T[i];}
 
private:
  void init() 
  {
    Z[0]=Z[1]=T[0]=T[1]=0.0;
    Nsol=0;
  }

  ROOT::Math::XYZTVector hita;
  ROOT::Math::XYZTVector hitb;

  //it can have 2 solutions
  double Z[2];
  double T[2];
  int Nsol;

};

void ZTcalc::compute()
{ 
  init();
  //T=A*Z+D
  double A=(hita.Z()-hitb.Z())/(hita.T()-hitb.T());
  double D=(hita.T()+hitb.T())/2-(hita.P2()-hitb.P2())/(hita.T()-hitb.T())/2;
  double B=D-hita.T();

  quadraticEquation q(1-A*A, -2*(hita.Z()+A*B), hita.P2()-B*B);

  Nsol=q.Nsolution();
  if (Nsol>=1)
    {
      Z[0]=q.solmoins();
      T[0]=A*Z[0]+D;
    }
  if (Nsol==2) 
    {
      Z[1]=q.solplus();
      T[1]=A*Z[1]+D;
    }

}
