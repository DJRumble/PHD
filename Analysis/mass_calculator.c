#include <stdio.h>
#include <math.h>


//Set Global constants

#define h 6.626068E-34
#define c_c 299792458
#define k 1.3806503E-23
#define psc 3.08567758E16
#define M_s 1.989E30



int main (void)
{

//deffine Parameters

  double m,M,S,D,A,B,nu;
  double lambda,omega,kappa,d,t,del;

  //deffine debugging parameters

  //double u,v,w,x,y,z;

  //set known parameters

  lambda = 850e-6; //m
  d = 130.0; //Pcs
  omega = 238.23; //arc secs ^2
  t = 10.0; //K
  kappa = 0.00191; //M^2/Kg
  del = 0.0016666666666; // CDEL2 from JCMT website degs

  //deffine input FLUX in Jy per beam

  S = 794.83E-26; 

  //find frequency - nu in Hz
 
  nu = c_c / lambda;
  printf("nu = %e\n",nu);
  
  //find distance - D in m

  D = d * psc;
  printf("D = %e\n",D);

  //find Black Body comp - B WHz^-1m^2

  //Start debugging//

  //x = ((2.0*h*pow(nu,3.0))/pow(c_c,2.0));
  //z = (h*nu)/(k*t);
  //y = (1.0/(exp(z)-1.0));
  //B = x*y;
  //printf("x = %e\n",x);
  //printf("z = %e\n",z);
  //printf("y = %e\n",y);

  //End debugging//

  B = (2.0*h*pow(nu,3.0))/(pow(c_c,2.0)*(1.0/(exp((h*nu)/(k*t))-1.0)));
  printf("blackbody = %e\n",B);

  //find beam correction - A (no unit)

  A = pow((del*3600.0),2.0)/omega;
  printf("correction = %e\n",A);

  // FIND MASS - M in solar masses

  //Start debugging//

  //u = (S*A*pow(D,2.0));
  //v = (B*kappa);
  //m = u/v;
  //printf("u = %e\n",u);
  //printf("v = %e\n",v);

  //End debugging//

  m = ((S*A*pow(D,2.0))/(B*kappa)); 
  M = m/M_s;

  printf("mass in kg = %e\n",m);
  printf("mass in solar masses = %e\n",M);
}
  
