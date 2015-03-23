#ifndef MAIN
#define MAIN

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "matrix.h"

//screen location and size constants
#define SX 512 //1536
#define SY 512
#define SCREENX 512
#define SCREENY 512
#define POSX 10
#define POSY 10

//camera variables
#define ASPECT_RATIO 1
#define VIEW_ANGLE 30
#define NEAR_PLANE 20
#define FAR_PLANE 300
//vectors
#define Ex 40.0
#define Ey 0.0
#define Ez 3.0
#define Gx 0.0
#define Gy 0.0
#define Gz 3.0
#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

//colour (RGB)
typedef struct{
	double r;
	double g;
	double b;
}Col;

//object 
typedef struct{
	char type;
	dmatrix_t trans;
	dmatrix_t invtrans;
	Col ambientCol;
	Col diffuseCol;
	Col specCol;
	double coeffs[3];
}Object;

//light source (and ambient)
typedef struct{
	dmatrix_t pos;
	Col ambient;
	Col col;
}Light;

//initialize objects
void readObjects();
// initialization camera
void initializeCamera();
//camera transformation functions
void generateCamMat();
void rotU(int);
void rotV(int);
void rotN(int);
//ray tracing 
void traceRays(Display *disp, Window win, int s);
void traceRay(Display *disp, Window win, int s, int i, int j);
int getMinT(dmatrix_t *e,dmatrix_t *d, Object *objectIntersect, dmatrix_t *pointIntersect);
int checkShadow(dmatrix_t point);
void thresholdLight(double *intensity);
//X11 display variables
Display *InitX(Display *disp, Window *win, int *s);
void SetCurrentColorX(Display *disp, GC *gc, unsigned int r, unsigned int g, unsigned int b);
void SetPixelX(Display *disp, Window win, int s, int i, int j);
void QuitX(Display *disp, Window win);
//vector and math functions
double vec3Dot(dmatrix_t *a, dmatrix_t *b);
double vec4Dot(dmatrix_t *a, dmatrix_t *b);
double vecMag(dmatrix_t *a, int elems);
double max(double a, double b);

#endif