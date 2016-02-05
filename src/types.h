#ifndef TYPES_H
#define TYPES_H

typedef float real;
typedef int t_scale;

//#define SSE 1

extern double const d2i;
#define qint(x,res)	*(int*)&(res=x+d2i)

#define HALFPI 1.57079632679490
#define PI 3.14159265358979
#define TWOPI 6.28318530717959
#define FOURPI 12.56637061435917


#endif
