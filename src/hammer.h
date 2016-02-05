#ifndef HAMMER_H
#define HAMMER_H
#include <stdio.h>

class Hammer {
public:
  Hammer(float Fs);
  void set(int upsample, float m, float K, float p, float Z, float alpha);
  ~Hammer();
  
  float getX();
  void strike(float v);
  float load(float vin, float vin1);
  bool isEscaped();
  
protected:
  float dt;
  float dti;
  float x;
  float v;
  float a;

  bool bEscaped;
  int escapeCount;
  int upsample;
  float mi;
  float K;
  float p;
  float Fs;
  float F;
  float upprev;
  float alpha;
  float Z2i;
};

#endif
