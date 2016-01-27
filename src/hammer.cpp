#include "hammer.h"
#include "types.h"
#include <math.h>
#include <stdio.h>

enum {
  S = 3
};

Hammer :: Hammer(float Fs)
{
  this->Fs = Fs;
  this->dt = 1.0/(Fs*S);
  this->dti = 1.0/dt;
  this->x = 0;
  this->a = 0.0;
  this->F = 0.0;
  this->upprev = 0;
}

void Hammer :: set(float m, float K, float p, float Z, float alpha) 
{
  this->p = p;
  this->K = K;
  this->mi = 1.0/m;
  this->alpha = alpha;
  this->Z2i = 1.0/(2.0*Z);
}

Hammer :: ~Hammer()
{
}

void Hammer :: strike(float v)
{
  this->x = 0;
  this->v = v;
}

float Hammer :: load(float vin) {
  static FILE *fp = fopen("hammer","w");
  for(int j=0;j<S;j++) {
    float up;
    up = (x>0)?pow((float)x,(float)p):0;
    float dupdt = (up-upprev)*dti;
    float v1;
		float x1;
		for(int k=0;k<S;k++) {
			F = K*(up+alpha*dupdt);
			if(F<=0) {
        F = 0;
      }       
			a = -F*mi;
			v1 = v + a * dt;
      x1 = x + (v1-(vin+F*Z2i)) * dt;
      
			float upnew = (x1>0)?pow((float)x1,(float)p):0;
			float dupdtnew = (upnew - upprev)*0.5*dti;
			float change = dupdtnew - dupdt;
			dupdt = dupdt + (float)0.5*change;
		} 
		v = v1;
		x = x1;		
		upprev = up;
  }
  
  //fprintf(fp,"%g\n",F);
  return F;
}

float Hammer :: getX()
{
  return x;
}
