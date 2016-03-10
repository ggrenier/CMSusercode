#ifndef ZvertexCalculator_HH
#define ZvertexCalculator_HH


#include <cmath>


//
// helper class to compute Z-vertex
//
struct hitData
{
  double x,y;
  double r;
  double z;
  double t;
  double localdirEta,localdirPhi,localdirTheta;
  int region;
  int station;
  int ring;
};



class ZvertexCalculator
{
public:
  ZvertexCalculator() : region(1), b(1), c(0), cosbeta(1), sinbeta(0), hit_time(0) {}


  inline ZvertexCalculator(hitData &hit) 
  {
    setHitData(hit);
  }

  void setHitData(hitData &hit);
  void setTimeShift(double timeshift);
  
  inline int computeSignOfCos()
  {
    if ((region==1 && b<c)||(region==-1 && b>c)) return 1;
    else return -1;
  }

  inline bool hasSolution() {return b>c*sinbeta;}
  inline double aOneSolution(){return sqrt(b*b-c*c*sinbeta*sinbeta)+c*cosbeta;}
  inline double aSecondSolution(){return c*cosbeta-sqrt(b*b-c*c*sinbeta*sinbeta);}
  inline bool hasTwoSolution() {return hasSolution() && b<c && cosbeta>0;}
  double aMinSolution();

private:
  //see http://en.wikipedia.org/wiki/Solution_of_triangles#Two_sides_and_non-included_angle_given_.28SSA.29
  //for notation
  int region;
  double b;
  double c;
  double cosbeta;
  double sinbeta;
  double hit_time;
};

#endif
