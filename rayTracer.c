#include "rayTracer.h"


//global camera vectors
dmatrix_t camU, camV, camN, camPos;
//camera matrix and inverse camera matrix
dmatrix_t cMat, ciMat;
//rotation matrix (for rotating and tilting camera)
dmatrix_t rotMat;

//global light source variable
Light light;

//Height and Width of screen
double H, W;

//Objects
Object obj[10];
int numObjects;

//pixelArray
double pixI[3];
double pixelArray[SX][SY][3];

//*** MAIN ***
//Saad Ahmed
//Nov 17, 2014
//main function
int main(void){
	//X11 variables
	Display *disp;
	Window win;
	int s;
	XEvent evt;
	KeySym key;
  	char text[255];
  	//for-loop counters
  	int i, j;

  	//initialize objects in the scene
  	readObjects();

  	//initialize camera matrix (Mv)
	initializeCamera();

  	//set up display and process user input
	disp = InitX(disp, &win, &s); 
	while (1){
		//get user events
	  	XNextEvent(disp, &evt);

	  	//fill background  on expose
	  	SetCurrentColorX(disp, &(DefaultGC(disp, s)), 0, 0, 0);
	    if (evt.type == Expose){
			for (i = 0; i < SCREENX; i+=1){
				for (j = 0; j < SCREENY; j+=1){
		          	SetPixelX(disp, win, s, i, j);
		    	}
		  	}
	    }

	    //handle user input
	    if(evt.type == KeyPress && XLookupString(&evt.xkey,text,255,&key,0)){
	    	switch(text[0]){
	    		//exit on ESC
		    	case 27: QuitX(disp,win);exit(0);
		    	//adjusting camera movements
		    	case 'w': camPos=*dmat_add(&camPos, &camV);camPos.m[4][1]=1;break;
		    	case 'a': camPos=*dmat_sub(&camPos, &camU);camPos.m[4][1]=1;break;
		    	case 's': camPos=*dmat_sub(&camPos, &camV);camPos.m[4][1]=1;break;
		    	case 'd': camPos=*dmat_add(&camPos, &camU);camPos.m[4][1]=1;break;
		    	//zoom in/out
		    	case 'e': camPos=*dmat_sub(&camPos, &camN);camPos.m[4][1]=1;break;
		    	case 'q': camPos=*dmat_add(&camPos, &camN);camPos.m[4][1]=1;break;
		    	//camera rotations and tilt
		    	case 'i': rotU(1);break;
		    	case 'j': rotV(1);break;
		    	case 'k': rotU(-1);break;
		    	case 'l': rotV(-1);break;
		    	case 'o': rotN(-1);break;
		    	case 'u': rotN(1);break;
	    	}

	    	traceRays(disp, win, s);
	    }
	}

	//exit
	QuitX(disp,win);
	return 0;
}

//*** readObjects ***
//Saad Ahmed
//Nov 17, 2014
//reads in all the objects from input and stores them 
void readObjects(){
	
	int i;	//for loop counter

	//calculate W and H
	H = NEAR_PLANE*tan((VIEW_ANGLE/2)*M_PI/180);
	W = H*ASPECT_RATIO;

	//define light position
	dmat_alloc(&(light.pos), 4, 1);
	scanf("%lf %lf %lf\n", &(light.pos.m[1][1]), &(light.pos.m[2][1]), &(light.pos.m[3][1]));
	light.pos.m[4][1]=1;
	//ambient light
	scanf("%lf %lf %lf\n", &(light.ambient.r), &(light.ambient.g), &(light.ambient.b));
	//point light colour intensity
	scanf("%lf %lf %lf\n", &(light.col.r), &(light.col.g), &(light.col.b));

	//other objects
	scanf("%d\n", &numObjects);
	for(i=0;i<numObjects; i++){
		//scan type
		scanf("%c\n", &(obj[i].type));
		//scan transformation matrix
		dmat_alloc(&(obj[i].trans), 4, 4);
		scanf("%lf %lf %lf %lf\n", &(obj[i].trans.m[1][1]), &(obj[i].trans.m[1][2]), &(obj[i].trans.m[1][3]), &(obj[i].trans.m[1][4]));
		scanf("%lf %lf %lf %lf\n", &(obj[i].trans.m[2][1]), &(obj[i].trans.m[2][2]), &(obj[i].trans.m[2][3]), &(obj[i].trans.m[2][4]));
		scanf("%lf %lf %lf %lf\n", &(obj[i].trans.m[3][1]), &(obj[i].trans.m[3][2]), &(obj[i].trans.m[3][3]), &(obj[i].trans.m[3][4]));
		scanf("%lf %lf %lf %lf\n", &(obj[i].trans.m[4][1]), &(obj[i].trans.m[4][2]), &(obj[i].trans.m[4][3]), &(obj[i].trans.m[4][4]));
		obj[i].invtrans=*dmat_inverse(&(obj[i].trans));
		//scan colours
		scanf("%lf %lf %lf\n", &(obj[i].ambientCol.r), &(obj[i].ambientCol.g), &(obj[i].ambientCol.b));
		scanf("%lf %lf %lf\n", &(obj[i].diffuseCol.r), &(obj[i].diffuseCol.g), &(obj[i].diffuseCol.b));
		scanf("%lf %lf %lf\n", &(obj[i].specCol.r), &(obj[i].specCol.g), &(obj[i].specCol.b));
		//scan coefficients
		scanf("%lf %lf %lf\n", &(obj[i].coeffs[0]), &(obj[i].coeffs[1]), &(obj[i].coeffs[2]));
	}
}

//*** initializeCamera ***
//Saad Ahmed
//Oct. 4, 2014
//function sets the values for vectors u, v, n, and e
void initializeCamera(){

	//calculate E, U, V, and N
	dmatrix_t E ; /* The center of projection for the camera */
    dmat_alloc(&E,4,1) ;
    E.m[1][1] = Ex ;
    E.m[2][1] = Ey ;
    E.m[3][1] = Ez ;
    E.m[4][1] = 1.0 ;
    dmatrix_t G ; /* Point gazed at by camera */
    dmat_alloc(&G,4,1) ;
    G.m[1][1] = Gx ;
    G.m[2][1] = Gy ;
    G.m[3][1] = Gz ;
    G.m[4][1] = 1.0 ;
    dmatrix_t N ;
    N = *dmat_normalize(dmat_sub(&E,&G)) ;
    N.l = 3 ;
    dmatrix_t UP ;
    dmat_alloc(&UP,4,1) ;
    UP.l = 3 ;
    UP.m[1][1] = UPx ;
    UP.m[2][1] = UPy ;
    UP.m[3][1] = UPz ;
    UP.m[4][1] = 1.0 ;
    dmatrix_t U ;
    U = *dmat_normalize(dcross_product(&UP,&N)) ;
    dmatrix_t V ;
    V = *dcross_product(&N,&U) ;

    //fill in camera's E, U, V and N
    dmat_alloc(&camPos, 4, 1);
	camPos.m[1][1]=E.m[1][1];camPos.m[2][1]=E.m[2][1];camPos.m[3][1]=E.m[3][1];camPos.m[4][1]=1;
	dmat_alloc(&camU, 4, 1);
	camU.m[1][1]=U.m[1][1];camU.m[2][1]=U.m[2][1];camU.m[3][1]=U.m[3][1];camU.m[4][1]=1;
	dmat_alloc(&camV, 4, 1);
	camV.m[1][1]=V.m[1][1];camV.m[2][1]=V.m[2][1];camV.m[3][1]=V.m[3][1];camV.m[4][1]=1;
	dmat_alloc(&camN, 4, 1);
	camN.m[1][1]=N.m[1][1];camN.m[2][1]=N.m[2][1];camN.m[3][1]=N.m[3][1];camN.m[4][1]=1;
	
	//final camera matrix and camera's rotation matrix
	dmat_alloc(&cMat,4,4);
	dmat_alloc(&rotMat,4,4);
}

//*** rotU ***
//Saad Ahmed
//Oct. 4, 2014
//function creates a rotation matrix about vector U
//and rotates vectors V and N with it
void rotU(int f){
	double SIN =sin(f*M_PI/180);
	double COS =cos(f*M_PI/180);
	double ux=camU.m[1][1], uy=camU.m[2][1], uz=camU.m[3][1];
	rotMat.m[1][1]=COS+ux*ux*(1-COS);rotMat.m[1][2]=ux*uy*(1-COS)-uz*SIN;rotMat.m[1][3]=ux*uz*(1-COS)+uy*SIN;rotMat.m[1][4]=0;
	rotMat.m[2][1]=uy*ux*(1-COS)+uz*SIN;rotMat.m[2][2]=COS+uy*uy*(1-COS);rotMat.m[2][3]=uy*uz*(1-COS)-ux*SIN;rotMat.m[2][4]=0;
	rotMat.m[3][1]=uz*ux*(1-COS)-uy*SIN;rotMat.m[3][2]=uz*uy*(1-COS)+ux*SIN;rotMat.m[3][3]=COS+uz*uz*(1-COS);rotMat.m[4][4]=0;
	rotMat.m[4][1]=0;rotMat.m[4][2]=0;rotMat.m[4][3]=0;rotMat.m[4][4]=1;
	camV = * dmat_mult(&rotMat, &camV);
	camN = * dmat_mult(&rotMat, &camN);
}

//*** rotV ***
//Saad Ahmed
//Oct. 4, 2014
//function creates a rotation matrix about vector U
//and rotates vectors U and N with it
void rotV(int f){
	double SIN =sin(f*M_PI/180);
	double COS =cos(f*M_PI/180);
	double ux=camV.m[1][1], uy=camV.m[2][1], uz=camV.m[3][1];
	rotMat.m[1][1]=COS+ux*ux*(1-COS);rotMat.m[1][2]=ux*uy*(1-COS)-uz*SIN;rotMat.m[1][3]=ux*uz*(1-COS)+uy*SIN;rotMat.m[1][4]=0;
	rotMat.m[2][1]=uy*ux*(1-COS)+uz*SIN;rotMat.m[2][2]=COS+uy*uy*(1-COS);rotMat.m[2][3]=uy*uz*(1-COS)-ux*SIN;rotMat.m[2][4]=0;
	rotMat.m[3][1]=uz*ux*(1-COS)-uy*SIN;rotMat.m[3][2]=uz*uy*(1-COS)+ux*SIN;rotMat.m[3][3]=COS+uz*uz*(1-COS);rotMat.m[4][4]=0;
	rotMat.m[4][1]=0;rotMat.m[4][2]=0;rotMat.m[4][3]=0;rotMat.m[4][4]=1;
	camU = * dmat_mult(&rotMat, &camU);
	camN = * dmat_mult(&rotMat, &camN);
}

//*** rotN ***
//Saad Ahmed
//Oct. 4, 2014
//function creates a rotation matrix about vector U
//and rotates vectors V and U with it
void rotN(int f){
	double SIN =sin(f*M_PI/180);
	double COS =cos(f*M_PI/180);
	double ux=camN.m[1][1], uy=camN.m[2][1], uz=camN.m[3][1];
	rotMat.m[1][1]=COS+ux*ux*(1-COS);rotMat.m[1][2]=ux*uy*(1-COS)-uz*SIN;rotMat.m[1][3]=ux*uz*(1-COS)+uy*SIN;rotMat.m[1][4]=0;
	rotMat.m[2][1]=uy*ux*(1-COS)+uz*SIN;rotMat.m[2][2]=COS+uy*uy*(1-COS);rotMat.m[2][3]=uy*uz*(1-COS)-ux*SIN;rotMat.m[2][4]=0;
	rotMat.m[3][1]=uz*ux*(1-COS)-uy*SIN;rotMat.m[3][2]=uz*uy*(1-COS)+ux*SIN;rotMat.m[3][3]=COS+uz*uz*(1-COS);rotMat.m[4][4]=0;
	rotMat.m[4][1]=0;rotMat.m[4][2]=0;rotMat.m[4][3]=0;rotMat.m[4][4]=1;
	camV = * dmat_mult(&rotMat, &camV);
	camU = * dmat_mult(&rotMat, &camU);
}

//*** generateCamMat ***
//Saad Ahmed
//Oct. 4, 2014
//function generates a camera transformation matrix from specified
//u, v, n, and e values
void generateCamMat(){
	cMat.m[1][1]=camU.m[1][1];cMat.m[1][2]=camU.m[2][1];cMat.m[1][3]=camU.m[3][1];
	cMat.m[2][1]=camV.m[1][1];cMat.m[2][2]=camV.m[2][1];cMat.m[2][3]=camV.m[3][1];
	cMat.m[3][1]=camN.m[1][1];cMat.m[3][2]=camN.m[2][1];cMat.m[3][3]=camN.m[3][1];
	cMat.m[4][1]=0;cMat.m[4][2]=0;cMat.m[4][3]=0;
	cMat.m[1][4]=-1*vec3Dot(&camU, &camPos);
	cMat.m[2][4]=-1*vec3Dot(&camV, &camPos);
	cMat.m[3][4]=-1*vec3Dot(&camN, &camPos);
	cMat.m[4][4]=1;
}

//*** thresholdLight ***
//Saad Ahmed
//Nov 17, 2014
//thresholds the light intensity so that values are between 0 and 1
void thresholdLight(double *intensity){
	if(intensity[0]<0) intensity[0]=0;
	if(intensity[0]>1) intensity[0]=1;
	if(intensity[1]<0) intensity[1]=0;
	if(intensity[1]>1) intensity[1]=1;
	if(intensity[2]<0) intensity[2]=0;
	if(intensity[2]>1) intensity[2]=1;
}

//*** traceRays ***
//Saad Ahmed
//Nov 14, 2014
//traces camera and shadow rays for each pixel on the screen
void traceRays(Display *disp, Window win, int s){
	int i, j;
	double currR, currG, currB;
	double mid = 1.0/9, out=1.0/9;

  	//trace a ray for each pixel and store intensity
  	for(i=0; i<SX; i++){
  		for(j=0; j<SY; j++){
			traceRay(disp, win, s, i, j);
  			pixelArray[i][j][0]=pixI[0];
  			pixelArray[i][j][1]=pixI[1];
  			pixelArray[i][j][2]=pixI[2];
  			thresholdLight(pixelArray[i][j]);
  		}
  		printf("%d\n", i);
  	}

  	//transfer onto screen (super resolution)
  	// for(i=0; i<SX; i+=3){
  	// 	for(j=0; j<SY; j+=3){
  	// 		currR=(pixelArray[i][j][0]+pixelArray[i+1][j][0]+pixelArray[i+2][j][0]+
  	// 				pixelArray[i][j+1][0]+pixelArray[i+1][j+1][0]+pixelArray[i+2][j+1][0]+
  	// 				pixelArray[i][j+2][0]+pixelArray[i+1][j+2][0]+pixelArray[i+2][j+2][0])/9.0;
  	// 		currG=(pixelArray[i][j][1]+pixelArray[i+1][j][1]+pixelArray[i+2][j][1]+
  	// 				pixelArray[i][j+1][1]+pixelArray[i+1][j+1][1]+pixelArray[i+2][j+1][1]+
  	// 				pixelArray[i][j+2][1]+pixelArray[i+1][j+2][1]+pixelArray[i+2][j+2][1])/9.0;
  	// 		currB=(pixelArray[i][j][2]+pixelArray[i+1][j][2]+pixelArray[i+2][j][2]+
  	// 				pixelArray[i][j+1][2]+pixelArray[i+1][j+1][2]+pixelArray[i+2][j+1][2]+
  	// 				pixelArray[i][j+2][2]+pixelArray[i+1][j+2][2]+pixelArray[i+2][j+2][2])/9.0;
  	// 		SetCurrentColorX(disp, &(DefaultGC(disp, s)), 
  	// 			(int)(currR*255), (int)(currG*255), (int)(currB*255));
  	// 		SetPixelX(disp, win, s, i/3, (SY-j)/3);
  	// 	}
  	// }
  	//transfer onto screen (non-super resolution)
  	for(i=0; i<SX; i++){
  		for(j=0; j<SY; j++){
  			currR=pixelArray[i][j][0];currG=pixelArray[i][j][1];currB=pixelArray[i][j][2];
  			SetCurrentColorX(disp, &(DefaultGC(disp, s)), 
  				(int)(currR*255), (int)(currG*255), (int)(currB*255));
  			SetPixelX(disp, win, s, i, SY-j);
  		}
  	}
}

//*** traceRay ***
//Saad Ahmed
//Nov 17, 2014
//given a column and row, this function traces rays for the pixel
//at that column and row
void traceRay(Display *disp, Window win, int s, int col, int row){
	int intersected;
	Object objectIntersect;
	double distanceToLight;
	double fallout;
	dmatrix_t pointIntersect, vecS, vecN, vecV, vecR, temp, vecReverseS;
	double diffDot, specDot;
	dmatrix_t *temp1, *temp2, *temp3, *temp4, *temp5;

	//calculate e and d for the ray, based on column and row
	dmatrix_t e = camPos;
	temp1=dmat_mult_constant(&camU, W*(2.0*col/SX-1));
	temp2=dmat_mult_constant(&camV, H*(2.0*row/SY-1));
	temp3=dmat_mult_constant(&camN, -1*NEAR_PLANE);
	temp4=dmat_add(temp3,temp1);
	dmatrix_t d = *dmat_add(temp4,temp2);
	d.m[4][1]=0;
	//free temps
	free_dmatrix(temp1->m, 1, temp1->l, 1, temp1->c);free_dmatrix(temp2->m, 1, temp2->l, 1, temp2->c);free(temp1);free(temp2);
	free_dmatrix(temp3->m, 1, temp3->l, 1, temp3->c);free_dmatrix(temp4->m, 1, temp4->l, 1, temp4->c);free(temp3);free(temp4);

	//calculate  point and object of intersection
	intersected=getMinT(&e, &d, &objectIntersect, &pointIntersect);

	//colour pixel black and exit if no intersection
	if(!intersected){
		free_dmatrix(d.m, 1, d.l, 1, d.c);
		pixI[0]=objectIntersect.coeffs[0]*light.ambient.r;
		pixI[1]=objectIntersect.coeffs[0]*light.ambient.g;
		pixI[2]=objectIntersect.coeffs[0]*light.ambient.b;
		return;
	}

	//check for shadows at this location
	if(checkShadow(pointIntersect)){
		free_dmatrix(d.m, 1, d.l, 1, d.c);
		free_dmatrix(pointIntersect.m, 1, pointIntersect.l, 1, pointIntersect.c);
		pixI[0]=objectIntersect.coeffs[0]*objectIntersect.ambientCol.r*light.ambient.r;
		pixI[1]=objectIntersect.coeffs[0]*objectIntersect.ambientCol.g*light.ambient.g;
		pixI[2]=objectIntersect.coeffs[0]*objectIntersect.ambientCol.b*light.ambient.b;
		return;
	}

	//get vector S
	temp =*dmat_sub(&(light.pos), &pointIntersect);
	distanceToLight=vecMag(&temp, 3);
	temp.m[4][1]=0;
	vecS=*dmat_normalize(&temp);
	free_dmatrix(temp.m, 1, temp.l, 1, temp.c);
	//get vector V
	dmat_alloc(&vecV, 4, 1);
	vecV.m[1][1]=camN.m[1][1];vecV.m[2][1]=camN.m[2][1];vecV.m[3][1]=camN.m[3][1];vecV.m[4][1]=0;
	//get vector N
	if(objectIntersect.type=='s') {vecN=*dmat_mult(&(objectIntersect.invtrans), &pointIntersect);vecN.m[4][1]=0;}
	else if(objectIntersect.type=='p') {dmat_alloc(&vecN,4,1); vecN.m[1][1]=0;vecN.m[1][1]=0;vecN.m[3][1]=1;vecN.m[4][1]=0;}
	else if(objectIntersect.type=='c') {vecN=*dmat_mult(&(objectIntersect.invtrans), &pointIntersect);vecN.m[3][1]=0;vecN.m[4][1]=0;}
	else if(objectIntersect.type=='o') {vecN=*dmat_mult(&(objectIntersect.invtrans), &pointIntersect);vecN.m[3][1]=1-vecN.m[3][1];vecN.m[4][1]=0;}
	//get vector R
	vecReverseS=*dmat_mult_constant(&vecS, -1);
	temp1=dmat_mult_constant(&vecN, 2.0*vec3Dot(&vecS, &vecN)/(vecMag(&vecN,3)*vecMag(&vecN,3)));
	temp =*dmat_add(&vecReverseS, temp1);
	temp.m[4][1]=0;
	vecR=*dmat_normalize(&temp);
	vecR.m[4][1]=0;
	//free temps
	free_dmatrix(vecReverseS.m, 1, vecReverseS.l, 1, vecReverseS.c);
	free_dmatrix(temp.m, 1, temp.l, 1, temp.c);free_dmatrix(temp1->m, 1, temp1->l, 1, temp1->c);free(temp1);

	//fallout for light distance
	fallout=(0.01+0.01*distanceToLight+ 0.01*distanceToLight*distanceToLight);
	// fallout=1;

	//get diffuse dot product
	diffDot=max(0, (vec3Dot(&vecS, &vecN)));
	//get specular dot product
	specDot=max(0, pow(vec3Dot(&vecR, &vecV), 15));	

	//calculate RGB values of pixel based on its specular, ambient, and diffuse components
	pixI[0]=objectIntersect.coeffs[0]*objectIntersect.ambientCol.r*light.ambient.r+
			objectIntersect.coeffs[1]/fallout*objectIntersect.diffuseCol.r*light.col.r*diffDot+
			objectIntersect.coeffs[2]/fallout*objectIntersect.specCol.r*light.col.r*specDot;
	pixI[1]=objectIntersect.coeffs[0]*objectIntersect.ambientCol.g*light.ambient.g+
			objectIntersect.coeffs[1]/fallout*objectIntersect.diffuseCol.g*light.col.g*diffDot+
			objectIntersect.coeffs[2]/fallout*objectIntersect.specCol.g*light.col.g*specDot;
	pixI[2]=objectIntersect.coeffs[0]*objectIntersect.ambientCol.b*light.ambient.b+
			objectIntersect.coeffs[1]/fallout*objectIntersect.diffuseCol.b*light.col.b*diffDot+
			objectIntersect.coeffs[2]/fallout*objectIntersect.specCol.b*light.col.b*specDot;

	//free matrices
	free_dmatrix(d.m, 1, d.l, 1, d.c);
	free_dmatrix(pointIntersect.m, 1, pointIntersect.l, 1, pointIntersect.c);
	free_dmatrix(vecS.m, 1, vecS.l, 1, vecS.c);free_dmatrix(vecV.m, 1, vecV.l, 1, vecV.c);
	free_dmatrix(vecN.m, 1, vecN.l, 1, vecN.c);free_dmatrix(vecR.m, 1, vecR.l, 1, vecR.c);
	//return pixel intensity (RGB)
	return;
}

//*** getMinT ***
//Saad Ahmed
//Nov 17, 2014
//given a ray (e and d), this function calculates the point and object
//of intersection of that ray
int getMinT(dmatrix_t *orige, dmatrix_t *origd, Object *objectIntersect, dmatrix_t *pointIntersect){
	double t0, t1;
	double minT=10000000;
	double a, b, c;
	int i;
	dmatrix_t *d, *e, *temp1, *temp2, *temp3;
	int intersected=0;
	dmat_alloc(pointIntersect, 4, 1);
	double zplane;

	//go through all objects
	for(i=0;i<numObjects;i++){
		//transform ray to conform to the object
		e=dmat_mult(&(obj[i].invtrans), orige);
		temp1=dmat_mult(&(obj[i].invtrans), origd); 
		d=dmat_normalize(temp1);

		//it's a sphere, find t0 and t1
		if(obj[i].type=='s'){
			//calculate a, b, and c for quadratic equation
			a = vecMag(d, 3)*vecMag(d, 3);
			b = vec3Dot(e, d);
			c = vecMag(e, 3)*vecMag(e, 3) -1;
			//get minimum of t0 and t1
			if(b*b-a*c>=0) {
				t0=-1.0*b/a-sqrt(b*b-a*c)/a;
				t1=-1.0*b/a+sqrt(b*b-a*c)/a;
				if(t0<minT && t0>=0) {intersected=1;minT=t0;*objectIntersect=obj[i];}
				if(t1<minT && t1>=0) {intersected=1;minT=t1;*objectIntersect=obj[i];}
			}	
		} 
		//if cylinder, find t0 and t1
		else if(obj[i].type=='c'){
			//calculate a, b, and c for quadratic equation
			a = d->m[1][1] * d->m[1][1] + d->m[2][1] * d->m[2][1];
			b = e->m[1][1] * d->m[1][1] + e->m[2][1] * d->m[2][1];
			c = e->m[1][1] * e->m[1][1] + e->m[2][1] * e->m[2][1]-1;
			//get minimum of t0 and t1
			if(b*b-a*c>0) {
				t0=-1.0*b/a-sqrt(b*b-a*c)/a;
				t1=-1.0*b/a+sqrt(b*b-a*c)/a;
				if(t0<minT && t0>=0) {intersected=1;minT=t0;*objectIntersect=obj[i];}
				if(t1<minT && t1>=0) {intersected=1;minT=t1;*objectIntersect=obj[i];}
			}	
		} 
		//if cone, get t0 and t1
		else if(obj[i].type=='o'){
			//calculate a, b, and c for quadratic equation
			a = d->m[1][1] * d->m[1][1] + d->m[2][1] * d->m[2][1] - 0.25*d->m[3][1] * d->m[3][1];
			b = e->m[1][1] * d->m[1][1] + e->m[2][1] * d->m[2][1] + 0.25*d->m[3][1]*(1-e->m[3][1]);
			c = e->m[1][1] * e->m[1][1] + e->m[2][1] * e->m[2][1] - 0.25*(1-e->m[3][1])*(1-e->m[3][1]);
			//get minimum of t0 and t1
			if(b*b-a*c>0) {
				t0=-1.0*b/a-sqrt(b*b-a*c)/a;
				t1=-1.0*b/a+sqrt(b*b-a*c)/a;
				if(t0<minT && t0>=0) {intersected=1;minT=t0;*objectIntersect=obj[i];}
				if(t1<minT && t1>=0) {intersected=1;minT=t1;*objectIntersect=obj[i];}
			}	
		} 
		//it's a plane, find t 
		else if(obj[i].type=='p'){
			t0= -1*(e->m[3][1]/d->m[3][1]);
			if(t0<minT && t0>0){intersected=1;minT=t0;*objectIntersect=obj[i];}
		}
		

		//calcalate point of intersection if object hit
		if(intersected){
			free_dmatrix(pointIntersect->m, 1, pointIntersect->l, 1, pointIntersect->c);//free(pointIntersect);
			temp2=dmat_mult_constant(d, minT);
			temp3 = dmat_add(e, temp2);
			temp3->m[4][1]=1.0;
			*pointIntersect = *dmat_mult(&(obj[i].trans) , temp3);
			free_dmatrix(temp2->m, 1, temp2->l, 1, temp2->c);free(temp2);
			free_dmatrix(temp3->m, 1, temp3->l, 1, temp3->c);free(temp3);
		}

		free_dmatrix(temp1->m, 1, temp1->l, 1, temp1->c);free(temp1);
		free_dmatrix(d->m, 1, d->l, 1, d->c);free(d);
		free_dmatrix(e->m, 1, e->l, 1, e->c);free(e);
	}
	
	
	//return whether ray intersected with an object
	return intersected;	
}

//*** checkShadow ***
//Saad Ahmed
//Nov 17, 2014
//returns 1 or 0 based on whether the given point contains a shadow
//i.e. if a ray from the point to the light source hits another object
int checkShadow(dmatrix_t point){
	double t0, t1;
	dmatrix_t *d, *e, *origd, *orige, *temp1;
	int shadowflag;
	double a, b, c, distance;
	int i;

	//find ray (e and d) that goes from point to light
	orige = &point;
	temp1=dmat_sub(&(light.pos), &point);
	distance=vecMag(temp1, 4);
	temp1->m[4][1]=0;
	origd=dmat_normalize(temp1);
	free_dmatrix(temp1->m,1,temp1->l,1,temp1->c);
	
	//go through all possible objects this ray can hit
	for(i=0;i<numObjects;i++){
		//set shadow flag to 0
		shadowflag=0;
		//transform ray to conform to object being processed
		e=dmat_mult(&(obj[i].invtrans), orige);
		d=dmat_mult(&(obj[i].invtrans), origd); 

		//if object is a sphere
		if(obj[i].type=='s'){
			//calculate a, b, and c for quadratic equation
			a = vecMag(d, 3)*vecMag(d, 3);
			b = vec3Dot(e, d);
			c = vecMag(e, 3)*vecMag(e, 3) -1;
			//check intersection
			if(b*b-a*c>=0) {
				t0=-1.0*b/a-sqrt(b*b-a*c)/a;
				t1=-1.0*b/a+sqrt(b*b-a*c)/a;
				if(t0>0.01 || t1>0.01) shadowflag=1;
			}
			
		} 
		//if object is a cylinder
		else if(obj[i].type=='c'){
			//calculate a, b, and c for quadratic equation
			a = d->m[1][1] * d->m[1][1] + d->m[2][1] * d->m[2][1];
			b = e->m[1][1] * d->m[1][1] + e->m[2][1] * d->m[2][1];
			c = e->m[1][1] * e->m[1][1] + e->m[2][1] * e->m[2][1]-1;
			//get minimum of t0 and t1
			if(b*b-a*c>0) {
				t0=-1.0*b/a-sqrt(b*b-a*c)/a;
				t1=-1.0*b/a+sqrt(b*b-a*c)/a;
				if( (distance>t0 && t0>0.01) || (distance>t1 && t1>0.01)) shadowflag=1;
			}	
		} 
		//if object is a cone
		else if(obj[i].type=='o'){
			//calculate a, b, and c for quadratic equation
			a = d->m[1][1] * d->m[1][1] + d->m[2][1] * d->m[2][1] - 0.25*d->m[3][1] * d->m[3][1];
			b = e->m[1][1] * d->m[1][1] + e->m[2][1] * d->m[2][1] + 0.25*d->m[3][1]*(1-e->m[3][1]);
			c = e->m[1][1] * e->m[1][1] + e->m[2][1] * e->m[2][1] - 0.25*(1-e->m[3][1])*(1-e->m[3][1]);
			//get minimum of t0 and t1
			if(b*b-a*c>0) {
				t0=-1.0*b/a-sqrt(b*b-a*c)/a;
				t1=-1.0*b/a+sqrt(b*b-a*c)/a;
				if( (distance>t0 && t0>0.01) || (distance>t1 && t1>0.01)) shadowflag=1;
			}	
		} 
		//if object is a plane
		else if(obj[i].type=='p'){
			t0= -1*(e->m[3][1]/d->m[3][1]);
			if(t0>0.01) shadowflag=1;
		}
		


		//return shadow if hit
		free_dmatrix(d->m,1,d->l,1,d->c);free(d);
		free_dmatrix(e->m,1,e->l,1,e->c);free(e);
		if(shadowflag){
			free_dmatrix(origd->m,1,origd->l,1,origd->c);free(origd);
			return 1;
		}
	}

	//if no object was hit, return 0 (no shadow)
	free_dmatrix(origd->m,1,origd->l,1,origd->c);free(origd);
	return 0;
}

//*** vec4Mag ***
//Saad Ahmed
//Nov 6, 2014
//computes magnitude of a vector
double vecMag(dmatrix_t *a, int elems){
	double sumSquares=0;
	int i;

	for(i=1;i<=elems;i++){
		sumSquares+=a->m[i][1]*a->m[i][1];
	}

	return sqrt(sumSquares);
}

//*** vec3Dot ***
//Saad Ahmed
//Oct. 4, 2014
//computes dot product of two 3-dimensional vectors
double vec3Dot(dmatrix_t *a, dmatrix_t *b){
	return a->m[1][1]*b->m[1][1]+a->m[2][1]*b->m[2][1]+
			a->m[3][1]*b->m[3][1];
}

//*** vec4Dot ***
//Saad Ahmed
//Oct. 4, 2014
//computes dot product of two 4-dimensional vectors
double vec4Dot(dmatrix_t *a, dmatrix_t *b){
	return a->m[1][1]*b->m[1][1]+a->m[2][1]*b->m[2][1]+
			a->m[3][1]*b->m[3][1]+a->m[4][1]*b->m[4][1];
}

//*** max ***
//Saad Ahmed
//Nov 17, 2014
//return the maximum of two doubles
double max(double a, double b){
	return a>b?a:b;
}

//************************************************************
//The following functions are created by Prof. Beauchemin
//************************************************************

Display *InitX(Display *disp, Window *win, int *s) {
  disp = XOpenDisplay(NULL);
  if(disp == NULL) {
    printf("Cannot open display\n");
    exit(1) ;
  }
  *s = DefaultScreen(disp);
  *win = XCreateSimpleWindow(disp, RootWindow(disp, *s), POSX, POSY, SCREENX, SCREENY, 1, BlackPixel(disp, *s), WhitePixel(disp, *s));
  Atom delWindow = XInternAtom(disp, "WM_DELETE_WINDOW", 0);
  XSetWMProtocols(disp, *win, &delWindow, 1);
  XSelectInput(disp, *win, ExposureMask | KeyPressMask | ButtonPressMask);
  XMapWindow(disp, *win);
  return(disp);
}

void SetCurrentColorX(Display *disp, GC *gc, unsigned int r, unsigned int g, unsigned int b) {
	XSetBackground(disp, *gc, 0 << 16 | 0 << 8 | 0);
  XSetForeground(disp, *gc, r << 16 | g << 8 | b);
}

void SetPixelX(Display *disp, Window win, int s, int i, int j) {
	if(!(i<0||i>SCREENX|j<0||j>SCREENY)){
		XDrawPoint(disp, win, DefaultGC(disp, s), i, j);
	}	
}

void QuitX(Display *disp, Window win) {
  XDestroyWindow(disp,win);
  XCloseDisplay(disp);
}
