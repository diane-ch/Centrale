#include "mex.h"
#include "math.h"


/* La fonction de calcul */
void MEX_SAFT(double c, int Fs, int Nt, int Nel, double *p_u, int Nx, int Nz, double *p_x, double *p_z, double *p_A, double *p_O )
{
    double r;
    double tau;
    int tau_round;
    int ind;
    int i;
    int ix ;
    int iz;

    for (ix=0; ix<Nx;ix++){
       
        for (iz=0; iz<Nz;iz++){

            for (i=0; i<Nel;i++){

                r = sqrt(p_z[iz]*p_z[iz] +(p_x[ix]-p_u[i])*(p_x[ix]-p_u[i]));
                tau = 2*r / c;
                tau_round = floor(tau*Fs+0.5);
                ind = tau_round;
                
                if (ind<Nt){  
                    p_O[ix*Nz+iz] = p_O[ix*Nz+iz] + p_A[i*Nt+ind];
                }
            }
            
        } 
    }
}

/* La fonction d'appel => NE RIEN CHANGER DANS CETTE FONCTION */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Déclarations variables C */
    double c;
	int Fs;
	int Nt;
	int Nel; 
	double *p_u;
	int Nx;
	int Nz;
	double *p_x;
	double *p_z;
	double *p_A;
	int NoutSize;
	double *p_O;
    
    /* Affectation des entrées */
    c = mxGetScalar(prhs[0]);
    Fs = mxGetScalar(prhs[1]);
    Nt = mxGetScalar(prhs[2]);
    Nel = mxGetScalar(prhs[3]);
    p_u = mxGetPr(prhs[4]);
    Nx = mxGetScalar(prhs[5]);
    Nz = mxGetScalar(prhs[6]);
    p_x = mxGetPr(prhs[7]);
    p_z = mxGetPr(prhs[8]);
    p_A = mxGetPr(prhs[9]);
    
    /* Affectation des sorties */
    NoutSize = Nx*Nz;
    plhs[0] = mxCreateNumericMatrix((mwSize)NoutSize,1,mxDOUBLE_CLASS,0);
    p_O = mxGetData(plhs[0]);
    
    /* Appel de la fonction de calcul */
    MEX_SAFT(c,Fs,Nt,Nel,p_u,Nx,Nz,p_x,p_z,p_A,p_O);
}
        




