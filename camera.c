/*            PURPOSE : An example of matrix inversion using matrix.h
 
 PREREQUISITES : matrix.h
 
 */

#include "matrix.h"

#define Ex 15.0
#define Ey 15.0
#define Ez 10.0

#define Gx 0.0
#define Gy 0.0
#define Gz 1.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

int main() {
    
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
    printf("\n N is: \n") ;
    write_dmatrix(&N) ;
    
    
    
    dmatrix_t UP ;
    
    dmat_alloc(&UP,4,1) ;
    UP.l = 3 ;
    
    UP.m[1][1] = UPx ;
    UP.m[2][1] = UPy ;
    UP.m[3][1] = UPz ;
    UP.m[4][1] = 1.0 ;
    
    dmatrix_t U ;
    
    U = *dmat_normalize(dcross_product(&UP,&N)) ;
    
    printf("\n U is: \n") ;
    
    write_dmatrix(dmat_normalize(&U)) ;
    
    
    dmatrix_t V ;
    
    V = *dcross_product(&N,&U) ;
    
    printf("\n V is:\n") ;
    
    write_dmatrix(&V) ;
    
    
    /* Build matrix M_v */
    
    dmatrix_t Mv ;
    
    dmat_alloc(&Mv,4,4) ;
    
    Mv.m[1][1] = U.m[1][1] ; Mv.m[1][2] = U.m[2][1] ; Mv.m[1][3] = U.m[3][1] ; Mv.m[1][4] = -1.0*(E.m[1][1]*U.m[1][1] + E.m[2][1]*U.m[2][1] + E.m[3][1]*U.m[3][1]) ;
    
    Mv.m[2][1] = V.m[1][1] ; Mv.m[2][2] = V.m[2][1] ; Mv.m[2][3] = V.m[3][1] ; Mv.m[2][4] = -1.0*(E.m[1][1]*V.m[1][1] + E.m[2][1]*V.m[2][1] + E.m[3][1]*V.m[3][1]) ;
    
    Mv.m[3][1] = N.m[1][1] ; Mv.m[3][2] = N.m[2][1] ; Mv.m[3][3] = N.m[3][1] ; Mv.m[3][4] = -1.0*(E.m[1][1]*N.m[1][1] + E.m[2][1]*N.m[2][1] + E.m[3][1]*N.m[3][1]) ;
    
    Mv.m[4][1] = 0.0 ; Mv.m[4][2] = 0.0 ; Mv.m[4][3] = 0.0 ; Mv.m[4][4] = 1.0 ;
    
    printf("\n Mv is:\n") ;
    
    write_dmatrix(&Mv) ;
    
    
    
    /* This is the gaze point */
    
    dmatrix_t Pw ;
    
    dmat_alloc(&Pw,4,1) ;
    
    Pw.m[1][1] = 0.0 ;
    Pw.m[2][1] = 0.0 ;
    Pw.m[3][1] = 1.0 ;
    Pw.m[4][1] = 1.0 ;
    
    printf("\n The gaze point in camera coordinates is:\n") ;
    
    write_dmatrix(dmat_mult(&Mv,&Pw));
    
    
    dmatrix_t Ew ;
    
    dmat_alloc(&Ew,4,1) ;
    
    Ew.m[1][1] = 15.0 ;
    Ew.m[2][1] = 15.0 ;
    Ew.m[3][1] = 10.0 ;
    Ew.m[4][1] = 1.0 ;
    
    printf("\n The center of projection in camera coordinates is:\n") ;
    
    write_dmatrix(dmat_mult(&Mv,&Ew)) ;
    
    
#define N 5.0
#define F 25.0
    
    /* Build matrix M_p */
    
    dmatrix_t Mp ;
    
    dmat_alloc(&Mp,4,4) ;
    dmat_identity(&Mp) ;
    
    float a = -1.0*(F + N)/(F - N) ;
    float b = -2.0*(F*N)/(F - N) ;
    
    Mp.m[1][1] = N ;
    Mp.m[2][2] = N ;
    Mp.m[3][3] = a ;
    Mp.m[3][4] = b ;
    Mp.m[4][3] = -1.0 ;
    Mp.m[4][4] = 0.0 ;
    
    printf("\n The gaze point after prespective transformation is:\n") ;
    
    dmatrix_t Pv ;
    
    Pv = *dmat_mult(&Mp,dmat_mult(&Mv,&Pw)) ;
    write_dmatrix(&Pv) ;
    
    
    /* perform projective transformation */
    
    Pv.m[1][1] /= Pv.m[4][1] ;
    Pv.m[2][1] /= Pv.m[4][1] ;
    Pv.m[3][1] /= Pv.m[4][1] ;
    Pv.m[4][1] /= Pv.m[4][1] ;
    
    printf("\n The gaze point after perspective projection is:\n") ;
    
    write_dmatrix(&Pv) ;
    
    /* Build matrices T_1 and S_1. We need view angle and aspect ratio to do this */
    
#define THETA 45.0
#define ASPECT 1.0
    
    /* Work out coordinates of near plane corners */
    
    float top = N*tan(M_PI/180.0*THETA/2.0) ;
    float right = ASPECT*top ;
    float bottom = -top ;
    float left = -right ;
    
    printf("\n(top,right) = (%f,%f)\n",top,right) ;
    printf("\n(top,left) = (%f,%f)\n",top,left) ;
    printf("\n(bottom,right) = (%f,%f)\n",bottom,right) ;
    printf("\n(bottom,left) = (%f,%f)\n",bottom,left) ;
    
    /* Put (top,right) in vector form and transform it into the cannonical view volume */

    dmatrix_t topright, topleft, bottomright, bottomleft ;
    
    dmat_alloc(&topright,4,1) ;
    //dmat_alloc(&topleft,4,1) ;
    //dmat_alloc(&bottomright,4,1) ;
    //dmat_alloc(&bottomleft,4,1) ;
    
    topright.m[1][1] = right ; /* X coordinate of top right */
    topright.m[2][1] = top ; /* Y coordinate of top right */
    topright.m[3][1] = -N ; /* top right is defined on the near plane */
    topright.m[4][1] = 1.0 ; /* homogeneous coordinates */
    
    printf("\n Before perspective transformation, top right corner of near plane is:\n") ;
    
    write_dmatrix(&topright) ;
    
    dmatrix_t T1 ;
    
    dmat_alloc(&T1,4,4) ;
    
    T1 = *dmat_identity(&T1) ;
    
    T1.m[1][4] = -(right + left)/2.0 ;
    T1.m[2][4] = -(top + bottom)/2.0 ;
    
    dmatrix_t S1 ;
    
    dmat_alloc(&S1,4,4) ;
    
    S1 = *dmat_identity(&S1) ;
    
    S1.m[1][1] = 2.0/(right - left) ;
    S1.m[2][2] = 2.0/(top - bottom) ;
    
    printf("\n After perspective transformation, the top right corner of the near plane is:\n") ;
    
    topright = *dmat_mult(&S1,dmat_mult(&T1,dmat_mult(&Mp,&topright))) ;
    
    write_dmatrix(&topright) ;
    
    
    
    printf("\n After perspective division, the top right corner of the near plane is:\n") ;
    
    topright.m[1][1] /= topright.m[4][1] ;
    topright.m[2][1] /= topright.m[4][1] ;
    topright.m[3][1] /= topright.m[4][1] ;
    topright.m[4][1] /= topright.m[4][1] ;
    
    write_dmatrix(&topright) ;
    
    /* Build matrices T2, S2, and W. For this we need the window size in pixels */
    
#define W  512
#define H  512
    
    dmatrix_t T2 ;
    dmatrix_t S2 ;
    dmatrix_t W2 ;
    
    dmat_alloc(&T2,4,4) ;
    dmat_alloc(&S2,4,4) ;
    dmat_alloc(&W2,4,4) ;
    
    T2 = *dmat_identity(&T2) ;
    S2 = *dmat_identity(&S2) ;
    W2 = *dmat_identity(&W2) ;
    
    T2.m[1][4] = 1.0 ;
    T2.m[2][4] = 1.0 ;
    
    S2.m[1][1] = W/2.0 ;
    S2.m[2][2] = H/2.0 ;
    
    W2.m[2][2] = -1.0 ;
    W2.m[2][4] = (double)H ;
    
    printf("\n Matrix T2 is:\n") ;
    write_dmatrix(&T2) ;
    
    printf("\n Matrix S2 is:\n") ;
    write_dmatrix(&S2) ;
    
    printf("\n Matrix W2 is:\n") ;
    write_dmatrix(&W2) ;
    
    printf("\n Top right corner of near plane in pixel coordinates is\n") ;
    
    write_dmatrix(dmat_mult(&W2,dmat_mult(&S2,dmat_mult(&T2,&topright)))) ;
    
}
