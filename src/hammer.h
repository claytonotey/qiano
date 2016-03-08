#ifndef HAMMER_H
#define HAMMER_H

class Hammer {
public:
  Hammer();
  void set(float Fs, float m, float K, float p, float Z, float alpha, int escapeDelay);
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
  int escapeDelay;
  int escapeCount;
  float mi;
  float K;
  float p;
  float F;
  float upprev;
  float alpha;
  float Z2i;
};

#endif
