// MATLAB MEX function particle-in-cell (PIC) simulation of the IEC

#define EPS0        8.85e-12f   
#define QE          1.6e-19f 
#define MYCLOCK     CLOCK_MONOTONIC
#define PI          3.1415926535f        
#define TWOPI       6.2831853071f   
#define PIOVERTWO   1.5707963268f
#define TWOOVERPI   0.6366197723f   
#define CUDATIMERS  1  
#define MAXPERIODS  2048
#include "mex.h"    
#include <stdio.h>
#include <stdint.h>	
#include <stdlib.h>	
#include <time.h>
#include <curand.h>   
#include <curand_kernel.h>  
#include "headers\hostFunctions.h"
#include "headers\gpuErrchkMacro.cuh"
#include "headers\IndexMapper.cuh"      // used in findDensity, findAcceleration, and PEKE
#include "headers\findDensity.cuh" 
#include "headers\performCollisions.cuh"
#include "headers\findPotential.cuh"
#include "headers\findElectricField.cuh"
#include "headers\findAcceleration.cuh"
#include "headers\moveParticles.cuh"
#include "headers\PEKE.cuh"             // contains two functions, one for finding the total kinetic energy, the other for finding the total potential energy

// Device pointers 
static float 
    *d_px,  *d_pr,           // particle positions
    *d_vx,  *d_vr,  *d_vt,   // particle velocities
    *d_axE, *d_arE,          // particle accelerations due to E-field
    *d_axB, *d_arB,          // particle accelerations due to B-field (not actual acceleration, just mag field times q/m)
    *d_Ex,  *d_Er,           // grid e-fields
    *d_Bx,  *d_Br,           // grid b-fields
    *d_rho,
    *d_phi,
    *d_phiUpdate,
    *d_phiVec,
    *d_KE, *d_PE, *d_E,
    *d_PM1, *d_PM2, *d_PM3, *d_PM4, *d_PM5, *d_PM6, *d_PM7, 
    *d_EX1, *d_EX2, *d_EX3, *d_EX4, *d_EX5, *d_ER1, *d_ER2, *d_ER3, 
    *d_CellVolumesPM, *d_bmaxPM,
    *d_bDirichlet, *d_zOffPM;

static bool *d_hitWall;
static bool *hitWall;
static int *pP;
static curandState *d_state;
static int *d_Ic;                   // particle linear indices (index of the closest cell node, needed for collisions)
static int *d_cellResidents, *d_cellOccupants, *d_cellOccupantsMax, *d_pP;
static float *d_fusionRate, *d_partNs, *d_partAs;
static float *timeTemp, *tempX, *tempR, *tempT;
static int allocate = 1;
static int init = 1;

static float meanx, meanxSign, lastMeanxSign, sumPeriods;
static int periodsCompleted;
static int itGlobal;
static int meanP;
static float meanPfloat = 0.0f;
static float lastMeanPfloat = 0.0f;
static int lastMeanP = 0;
static float *E;
void cleanup(void) {
    mexPrintf("MEX-file is terminating, destroying array\n");
    cudaFree(d_px);
    cudaFree(d_pr);
    cudaFree(d_vx);
    cudaFree(d_vr);
    cudaFree(d_vt);
    cudaFree(d_axE);
    cudaFree(d_arE);
    cudaFree(d_axB);
    cudaFree(d_arB);
    cudaFree(d_Ex);
    cudaFree(d_Er);
    cudaFree(d_Bx);
    cudaFree(d_Br);
    cudaFree(d_rho);
    cudaFree(d_phi);
    cudaFree(d_phiUpdate);
    cudaFree(d_CellVolumesPM);
    cudaFree(d_bmaxPM);
    cudaFree(d_PM1);
    cudaFree(d_PM2);
    cudaFree(d_PM3);
    cudaFree(d_PM4);
    cudaFree(d_PM5);
    cudaFree(d_PM6);
    cudaFree(d_PM7);
    cudaFree(d_phiVec);
    cudaFree(d_EX1); 
    cudaFree(d_EX2);
    cudaFree(d_EX3);
    cudaFree(d_EX4);
    cudaFree(d_EX5);
    cudaFree(d_ER1);
    cudaFree(d_ER2);
    cudaFree(d_ER3);
    cudaFree(d_bDirichlet);
    cudaFree(d_zOffPM);
    cudaFree(d_KE);
    cudaFree(d_PE);
    cudaFree(d_E);
    cudaFree(d_hitWall);
    cudaFree(d_Ic); 
    cudaFree(d_cellResidents);
    cudaFree(d_partNs);
    cudaFree(d_partAs);
    cudaFree(d_cellOccupants);
    cudaFree(d_cellOccupantsMax);
    cudaFree(d_fusionRate);
    cudaFree(d_pP);
    free(hitWall);
    free(pP);
    free(tempX);
    free(tempR);
    free(tempT);
    mexPrintf("Freed CUDA Array.\n");
}

void kernelFunction(
        int load,
        float *px, float *pr, float *vx, float *vr, float *vt, float *Ex, float *Er, float *Bx, float *Br, float *rho, float *phi, float *phiVec,
        float *KE, float *PE,
        const float dx, const float dr, const int nx, const int nr, const int nEM, const int nPM,
        const int nxBound, const int nrBound, const int nInside, const float xBound, const float zOff,
        const float c, const float C, const float drSlope,
        const float dt, const int nt, const int numPeriods, int *np,
        const float creationRate, float *remainder, const float beamPoint, const float creationRadius, const float fusionVelocity,
        const float LxZero, const float LrZero, const float tanangle, const float q, const float m, const float qom,
        const size_t npBuffer,
        const float *PM1, const float *PM2, const float *PM3, const float *PM4, const float *PM5, const float *PM6, const float *PM7,
        const int PM_OFF_1A, const int PM_OFF_1B, const int PM_OFF_1T, const int PM_OFF_2A, const int PM_OFF_2B, const int PM_OFF_2T, const int PM_OFF_6A, const int PM_OFF_6B, const int PM_OFF_6T, const int PM_OFF_7A, const int PM_OFF_7B, const int PM_OFF_7T,
        float *EX1, float *EX2, float *EX3, float *EX4, float *EX5, float *ER1, float *ER2, float *ER3,
        const int EX_OFF_1A, const int EX_OFF_1B, const int EX_OFF_1T, const int EX_OFF_5A, const int EX_OFF_5B, const int EX_OFF_5T, const int ER_OFF_1A, const int ER_OFF_1B, const int ER_OFF_1T, const int ER_OFF_3A, const int ER_OFF_3B, const int ER_OFF_3T,
        const float *CellVolumesPM, const float *bmaxPM, const float a_coeff, const float macroWeight, const float W, const int itPM, const float L, const float *bDirichlet, const float *zOffPM,
        const int maxParticlesPerCell, int *cellResidents, float *partNs, float *partAs, float *fusionRate, int *cellOccupants, int extSeed, const bool FixEnergy, const bool loadPrevMode, const bool energyBump, float *Timers, float *energyLossTotal, float *fusionRateTotal,
        float *periods)
{
    
//     printf("Inside function\n");
    cudaSetDevice(0);
    const int blocks = 64;
    const int threads = 48;
    float Einitial, Etotal, KEtotal;
    
    if (allocate){
        // if memory has not been allocated, allocate it
        printf("allocating GPU memory\n");
        cudaMalloc((void**)&d_px,                npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_pr,                npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_vx,                npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_vr,                npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_vt,                npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_axE,               npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_arE,               npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_axB,               npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_arB,               npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_Ic,                npBuffer*               sizeof(int  ));
        cudaMalloc((void**)&d_Ex,                nEM*                    sizeof(float));
        cudaMalloc((void**)&d_Er,                nEM*                    sizeof(float));
        cudaMalloc((void**)&d_Bx,                nEM*                    sizeof(float));
        cudaMalloc((void**)&d_Br,                nEM*                    sizeof(float));
        cudaMalloc((void**)&d_rho,               nPM*                    sizeof(float));
        cudaMalloc((void**)&d_phi,               nPM*                    sizeof(float));
        cudaMalloc((void**)&d_phiUpdate,         nPM*                    sizeof(float));
        cudaMalloc((void**)&d_CellVolumesPM,     nPM*                    sizeof(float));
        cudaMalloc((void**)&d_bmaxPM,            nPM*                    sizeof(float));
        cudaMalloc((void**)&d_PM1,               nPM*                    sizeof(float));
        cudaMalloc((void**)&d_PM2,               nPM*                    sizeof(float));
        cudaMalloc((void**)&d_PM3,               nPM*                    sizeof(float));
        cudaMalloc((void**)&d_PM4,               nPM*                    sizeof(float));
        cudaMalloc((void**)&d_PM5,               nPM*                    sizeof(float));
        cudaMalloc((void**)&d_PM6,               nPM*                    sizeof(float));
        cudaMalloc((void**)&d_PM7,               nPM*                    sizeof(float));
        cudaMalloc((void**)&d_phiVec,            nEM*                    sizeof(float));
        cudaMalloc((void**)&d_EX1,               nEM*                    sizeof(float));
        cudaMalloc((void**)&d_EX2,               nEM*                    sizeof(float));
        cudaMalloc((void**)&d_EX3,               nEM*                    sizeof(float));
        cudaMalloc((void**)&d_EX4,               nEM*                    sizeof(float));
        cudaMalloc((void**)&d_EX5,               nEM*                    sizeof(float));
        cudaMalloc((void**)&d_ER1,               nEM*                    sizeof(float));
        cudaMalloc((void**)&d_ER2,               nEM*                    sizeof(float));
        cudaMalloc((void**)&d_ER3,               nEM*                    sizeof(float));
        cudaMalloc((void**)&d_bDirichlet,        nPM*                    sizeof(float));
        cudaMalloc((void**)&d_zOffPM,            nPM*                    sizeof(float));
        cudaMalloc((void**)&d_KE,                npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_PE,                npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_E,                 npBuffer*               sizeof(float));
        cudaMalloc((void**)&d_hitWall,           npBuffer*               sizeof(bool ));
        cudaMalloc((void**)&d_state,             nPM*                    sizeof(curandState));
        cudaMalloc((void**)&d_cellResidents,     maxParticlesPerCell*nPM*sizeof(int  ));
        cudaMalloc((void**)&d_partNs,            maxParticlesPerCell*nPM*sizeof(float));
        cudaMalloc((void**)&d_partAs,            maxParticlesPerCell*nPM*sizeof(float));
        cudaMalloc((void**)&d_fusionRate,        nPM*                    sizeof(float));
        cudaMalloc((void**)&d_cellOccupants,     nPM*                    sizeof(int));
        cudaMalloc((void**)&d_cellOccupantsMax,  nPM*                    sizeof(int));
        cudaMalloc((void**)&d_pP,                npBuffer*               sizeof(int));
        hitWall  =  (bool  *) malloc (npBuffer*sizeof (bool ));
        pP       =  (int   *) malloc (npBuffer*sizeof (int  ));
        timeTemp =  (float *) malloc (npBuffer*sizeof (float));
        tempX    =  (float *) malloc (npBuffer*sizeof (float));
        tempR    =  (float *) malloc (npBuffer*sizeof (float));
        tempT    =  (float *) malloc (npBuffer*sizeof (float));
        E        =  (float *) malloc (npBuffer*sizeof (float));
        allocate = 0;
        timeTemp[0] = 0.0f;
//         if(!loadPrevMode){
//             printf("calculating initial temperature\n");
//             FindTemperature(tempX,tempR,tempT,vx,vr,vt,0,npBuffer, fusionVelocity);
//         }
    }
     

    
    if(load){ // copy input vectors from host memory to GPU buffers
        printf("copying data from CPU to GPU\n");
        cudaMemcpy(d_px,              px,             npBuffer*               sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_pr,              pr,             npBuffer*               sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_vx,              vx,             npBuffer*               sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_vr,              vr,             npBuffer*               sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_vt,              vt,             npBuffer*               sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_Ex,              Ex,             nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_Er,              Er,             nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_Bx,              Bx,             nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_Br,              Br,             nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_PM1,             PM1,            nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_PM2,             PM2,            nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_PM3,             PM3,            nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_PM4,             PM4,            nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_PM5,             PM5,            nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_PM6,             PM6,            nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_PM7,             PM7,            nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_EX1,             EX1,            nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_EX2,             EX2,            nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_EX3,             EX3,            nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_EX4,             EX4,            nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_EX5,             EX5,            nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_ER1,             ER1,            nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_ER2,             ER2,            nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_ER3,             ER3,            nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_phi,             phi,            nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_phiVec,          phiVec,         nEM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_CellVolumesPM,   CellVolumesPM,  nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_bmaxPM,          bmaxPM,         nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_bDirichlet,      bDirichlet,     nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_zOffPM,          zOffPM,         nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_cellResidents,   cellResidents,  maxParticlesPerCell*nPM*sizeof(int  ), cudaMemcpyHostToDevice);
        cudaMemcpy(d_partNs,          partNs,         maxParticlesPerCell*nPM*sizeof(int  ), cudaMemcpyHostToDevice);
        cudaMemcpy(d_partAs,          partAs,         maxParticlesPerCell*nPM*sizeof(int  ), cudaMemcpyHostToDevice);
        cudaMemcpy(d_fusionRate,      fusionRate,     nPM*                    sizeof(float), cudaMemcpyHostToDevice);
        setupRandState<<<blocks,threads>>>(extSeed, nPM, d_state);

        cudaMemcpy(d_cellOccupants, cellOccupants, nPM*sizeof(int), cudaMemcpyHostToDevice);
        gpuErrchk(cudaMemset(d_pP, 0, npBuffer*sizeof(int))); // set max recorded number of each cell to zero

        periodsCompleted = 0;
        sumPeriods = 0;
        itGlobal = 0;
        mexAtExit(cleanup);

//         if(loadPrevMode){
            printf("Calculating Temp from CPUload\n");
            FindTemperature(tempX,tempR,tempT,vx,vr,vt,0,npBuffer, fusionVelocity);
//         }
    }
 
//     FindTemperature(tempX,tempR,tempT,vx,vr,vt,0,npBuffer, fusionVelocity);

    
    #if CUDATIMERS
    cudaEvent_t
            startCreate, stopCreate,
            startDensity, stopDensity,
            startPotential, stopPotential,
            startEfield, stopEfield,
            startAcceleration, stopAcceleration,
            startCollisions, stopCollisions,
            startMove, stopMove,
            startBoolTransfer, stopBoolTransfer,
            startKill, stopKill;
    
    cudaEventCreate(&startCreate);
    cudaEventCreate(&stopCreate);
    float timeCreate = 0;

    cudaEventCreate(&startDensity);
    cudaEventCreate(&stopDensity);
    float timeDensity = 0;

    cudaEventCreate(&startPotential);
    cudaEventCreate(&stopPotential);
    float timePotential = 0;
    
    cudaEventCreate(&startEfield);
    cudaEventCreate(&stopEfield);
    float timeEfield = 0;
    
    cudaEventCreate(&startAcceleration);
    cudaEventCreate(&stopAcceleration);
    float timeAcceleration = 0;
    
    cudaEventCreate(&startCollisions);
    cudaEventCreate(&stopCollisions);
    float timeCollisions = 0;
    
    cudaEventCreate(&startMove);
    cudaEventCreate(&stopMove);
    float timeMove = 0;
    
    cudaEventCreate(&startBoolTransfer);
    cudaEventCreate(&stopBoolTransfer);
    float timeBoolTransfer = 0;
    
    cudaEventCreate(&startKill);
    cudaEventCreate(&stopKill);
    float timeKill = 0;
    
    float tempTime = 0; 
    
    #endif
//     printf("2Bx[0] = %e\n",Bx[0]);

    clock_t startLoop = clock();
    clock_t startCreateCPU;
    float timeCreateCPU = 0;
    gpuErrchk(cudaMemset(d_cellOccupantsMax, 0, nPM*sizeof(int))); // set max recorded number of each cell to zero
    printf("starting loop with %i timesteps\n",nt);
    //for(int it = 0; it < nt; it++){
    int it = 0;
    while(periodsCompleted < numPeriods && it < nt && np[0]>0){

        // CALCULATE DENSITY
        #if CUDATIMERS  
        cudaEventRecord(startDensity);     
        #endif   
        gpuErrchk(cudaMemset(d_rho, 0, nPM*sizeof(float))); // set density to zero
        gpuErrchk(cudaMemset(d_cellOccupants, 0, nPM*sizeof(int))); // set number of particles in each cell to zero
        cudaDeviceSynchronize();     
        findDensity<<<blocks, threads>>>(  
                d_px, d_pr, 
                d_Ic,
                d_rho,
                dx, dr,
                nx-1, nr-1,
                nxBound-1, nrBound-1, xBound,
                c, C, drSlope, macroWeight,
                maxParticlesPerCell, d_cellResidents, d_cellOccupants,
                (*np));
        
        gpuErrchk(cudaPeekAtLastError()); gpuErrchk(cudaDeviceSynchronize());
        #if CUDATIMERS
        cudaEventRecord(stopDensity);
        cudaEventSynchronize(stopDensity);
        cudaEventElapsedTime(&tempTime, startDensity, stopDensity);
        timeDensity += tempTime;
        #endif
        
        // CALCULATE POTENTIAL
        #if CUDATIMERS
        cudaEventRecord(startPotential);
        #endif 
        for(int itP = 0; itP < itPM; itP++){ 
            findPotential<<<blocks, threads>>>(
                    d_rho, d_phi, d_phiUpdate, d_CellVolumesPM, d_bDirichlet, d_zOffPM,
                    d_PM1, d_PM2,d_PM3,d_PM4,d_PM5,d_PM6,d_PM7,
                    PM_OFF_1A,PM_OFF_1B,PM_OFF_1T,PM_OFF_2A,PM_OFF_2B,PM_OFF_2T,PM_OFF_6A,PM_OFF_6B,PM_OFF_6T,PM_OFF_7A,PM_OFF_7B,PM_OFF_7T,
                    nPM, W, L, q);
            findPotential<<<blocks, threads>>>(
                    d_rho, d_phiUpdate, d_phi, d_CellVolumesPM, d_bDirichlet, d_zOffPM,
                    d_PM1,d_PM2,d_PM3,d_PM4,d_PM5,d_PM6,d_PM7,
                    PM_OFF_1A,PM_OFF_1B,PM_OFF_1T,PM_OFF_2A,PM_OFF_2B,PM_OFF_2T,PM_OFF_6A,PM_OFF_6B,PM_OFF_6T,PM_OFF_7A,PM_OFF_7B,PM_OFF_7T,
                    nPM, W, L, q);
            gpuErrchk(cudaPeekAtLastError()); gpuErrchk(cudaDeviceSynchronize());
        }  
        // Copy unknown potentials from d_phi to d_phiVec - in the bulk region
        for (int ir = 0; ir < nrBound-1; ir++)
            gpuErrchk(cudaMemcpy(d_phiVec+ir*nx, d_phi+ir*(nx-1), (nx-1)*sizeof(float),cudaMemcpyDeviceToDevice));
        // Copy unknown potentials from d_phi to d_phiVec - first row above, where EM extends much farther than PM
        gpuErrchk(cudaMemcpy(d_phiVec+(nrBound-1)*nx, d_phi+(nrBound-1)*(nx-1), (nxBound-1)*sizeof(float),cudaMemcpyDeviceToDevice));
        // Copy unknown potentials from d_phi to d_phiVec - extra fusion region
        for (int ir = 0; ir < nr-nrBound-1; ir++)
            gpuErrchk(cudaMemcpy(d_phiVec+nx*nrBound+ir*nxBound, d_phi+(nrBound-1)*(nx-1)+(nxBound-1)+ir*(nxBound-1), (nxBound-1)*sizeof(float),cudaMemcpyDeviceToDevice));
        
        #if CUDATIMERS
        cudaEventRecord(stopPotential);
        cudaEventSynchronize(stopPotential);
        cudaEventElapsedTime(&tempTime, startPotential, stopPotential);
        timePotential += tempTime;
        #endif 

        // If on first iteration, calculate energy
        if (it == 1){
            findTotalEnergy<<<blocks, threads>>>(d_E,d_px,d_pr,d_vx,d_vr,d_vt,d_phi,dx,dr,nx,nr,nxBound,nrBound,nInside,xBound,c,C,drSlope,(*np),q,m,macroWeight);
            cudaMemcpy(E, d_E, (*np)*sizeof(float), cudaMemcpyDeviceToHost);
            Einitial = 0;
            for(int p = 0; p < (*np); p++){
                Einitial += E[p];
            }
        }
        
        // CALCULATE ELECTRIC FIELD
        #if CUDATIMERS
        cudaEventRecord(startEfield);
        #endif
        findElectricField<<<blocks, threads>>>( 
                d_phiVec, d_Ex, d_Er,
                d_EX1, d_EX2, d_EX3, d_EX4, d_EX5, d_ER1, d_ER2, d_ER3,
                EX_OFF_1A, EX_OFF_1B, EX_OFF_1T, EX_OFF_5A, EX_OFF_5B, EX_OFF_5T,
                ER_OFF_1A, ER_OFF_1B, ER_OFF_1T, ER_OFF_3A, ER_OFF_3B, ER_OFF_3T,
                nEM, L);
        gpuErrchk(cudaPeekAtLastError()); gpuErrchk(cudaDeviceSynchronize());
        #if CUDATIMERS
        cudaEventRecord(stopEfield);
        cudaEventSynchronize(stopEfield);
        cudaEventElapsedTime(&tempTime, startEfield, stopEfield);
        timeEfield += tempTime;
        #endif
        
        // INTERPOLATE ACCELERATION TO PARTICLES
        #if CUDATIMERS
        cudaEventRecord(startAcceleration);
        #endif

        findAcceleration<<<blocks, threads>>>(
                d_px, d_pr,
                d_axE, d_arE,
                d_axB, d_arB,
                d_Ex, d_Er,
                d_Bx, d_Br,
                dx, dr,
                nx, nr, 
                nxBound, nrBound, nInside, xBound,
                c, C, drSlope, qom, 
                (*np)); 
        gpuErrchk(cudaPeekAtLastError()); gpuErrchk(cudaDeviceSynchronize());

        #if CUDATIMERS
        cudaEventRecord(stopAcceleration);
        cudaEventSynchronize(stopAcceleration);
        cudaEventElapsedTime(&tempTime, startAcceleration, stopAcceleration);
        timeAcceleration += tempTime;
        #endif

        // PERFORM COLLISIONS
        #if CUDATIMERS
        cudaEventRecord(startCollisions);
        #endif
        
        // Then perform the collisions
        performCollisions<<<blocks, threads>>>( 
                d_vx, d_vr, d_vt,
                d_Ic,
                d_rho, d_bmaxPM, d_CellVolumesPM,
                maxParticlesPerCell, d_cellResidents, d_cellOccupants, d_cellOccupantsMax,
                d_fusionRate, d_partNs, d_partAs, 
                d_state,  
                PI*dt, a_coeff,
                nPM, L);
        

        #if CUDATIMERS
        cudaEventRecord(stopCollisions);
        cudaEventSynchronize(stopCollisions);
        cudaEventElapsedTime(&tempTime, startCollisions, stopCollisions);
        timeCollisions += tempTime;
        #endif
// 
        // -----------------------     
        // *** MOVE PARTICLES ***
        // -----------------------
        #if CUDATIMERS
        cudaEventRecord(startMove);
        #endif
        moveParticles<<<blocks, threads>>>(
                d_px, d_pr,
                d_vx, d_vr, d_vt,
                d_axE, d_arE,
                d_axB, d_arB,
                d_hitWall, d_pP, LxZero, LrZero, tanangle,
                xBound,
                (*np), dt);  
        gpuErrchk(cudaPeekAtLastError()); gpuErrchk(cudaDeviceSynchronize());

        #if CUDATIMERS
        cudaEventRecord(stopMove);
        cudaEventSynchronize(stopMove);
        cudaEventElapsedTime(&tempTime, startMove, stopMove);
        timeMove += tempTime;
        #endif 


        // TRANSFER BOOLEAN KILL LIST
        #if CUDATIMERS 
        cudaEventRecord(startBoolTransfer);
        #endif 
        cudaMemcpy(hitWall, d_hitWall, npBuffer*sizeof(bool), cudaMemcpyDeviceToHost);
        #if CUDATIMERS
        cudaEventRecord(stopBoolTransfer);
        cudaEventSynchronize(stopBoolTransfer);
        cudaEventElapsedTime(&tempTime, startBoolTransfer, stopBoolTransfer);
        timeBoolTransfer += tempTime;
        #endif
        
            findKineticEnergy<<<blocks, threads>>>(
                d_vx, d_vr, d_vt,
                d_KE,
                m, 
                (*np),
                macroWeight);
            
        float energyLossThisTimeStep = 0;
        for(int i = 0; i < (*np); i++){
            if(hitWall[i]){
                float particleEnergy;
                cudaMemcpy(&particleEnergy, d_KE+i, sizeof(float), cudaMemcpyDeviceToHost);
                energyLossThisTimeStep += particleEnergy;
            }
        }
                
        // KILL PARTICLES
        #if CUDATIMERS
        cudaEventRecord(startKill);
        #endif
        int npKill = 0;
        for (int p = 0; p < (*np); p++)
            if (hitWall[p])
                npKill++;
        np[1] += npKill;
        int npNew = (*np) - npKill;
        
        int pp = npNew-1;
        while (pp >= 0){
            if(hitWall[pp]){
                bool notFound = true;
                int ppp = (*np)-1; 
                while(notFound){
                    if(!(hitWall[ppp])){
                        cudaMemcpy(d_px+pp, d_px+ppp, sizeof(float), cudaMemcpyDeviceToDevice);
                        cudaMemcpy(d_pr+pp, d_pr+ppp, sizeof(float), cudaMemcpyDeviceToDevice);
                        cudaMemcpy(d_vx+pp, d_vx+ppp, sizeof(float), cudaMemcpyDeviceToDevice);
                        cudaMemcpy(d_vr+pp, d_vr+ppp, sizeof(float), cudaMemcpyDeviceToDevice);
                        cudaMemcpy(d_vt+pp, d_vt+ppp, sizeof(float), cudaMemcpyDeviceToDevice);
                        cudaMemcpy(d_pP+pp, d_pP+ppp, sizeof(int  ), cudaMemcpyDeviceToDevice);
                        notFound = false;
                        (*np)--;
                    }else{
                        ppp--;
                    }
                }
            }
            pp--;
        }
        (*np) = npNew;
        #if CUDATIMERS
        cudaEventRecord(stopKill);
        cudaEventSynchronize(stopKill);
        cudaEventElapsedTime(&tempTime, startKill, stopKill);
        timeKill += tempTime;
        #endif
        cudaMemcpy(fusionRate, d_fusionRate, nPM*sizeof(float), cudaMemcpyDeviceToHost);
//         for (int gah = 0; gah < nPM; gah++){
//             printf("fusion = %e\n",fusionRate[gah]) ;
//         }
        float fusionRateThisTimeStep = 0;
        for (int cc = 0; cc < nPM; cc++){
           fusionRateThisTimeStep+=fusionRate[cc];
        }
        fusionRateTotal[it] = fusionRateThisTimeStep;
        energyLossTotal[it] = energyLossThisTimeStep; 

        
        
        cudaMemcpy(pP, d_pP, (*np)*sizeof(int), cudaMemcpyDeviceToHost);
        meanP = 0;
        for(int p = 0; p < (*np); p++){
            meanP += pP[p];
        }
        meanPfloat = ((float)meanP)/((float)((*np)));
        float floorMeanPfloat = floor(meanPfloat);
        meanP = (int)(floorMeanPfloat);
        if (meanP > lastMeanP){
            printf("meanP = %i, lastMeanP = %i, meanPfloat = %f, lastMeanPfloat = %f\n",meanP, lastMeanP, meanPfloat, lastMeanPfloat);
            float dt_frac_backtrack = (lastMeanPfloat-floorMeanPfloat)/(meanPfloat-lastMeanPfloat);  // negative number that tells us how much of a timestep we overshot the period change by
            periods[periodsCompleted] = ((float)(itGlobal) + dt_frac_backtrack)*dt - sumPeriods;
            printf("completed period %i of %i in %e, dtFrac = %f, sumP = %e\n",periodsCompleted+1, numPeriods, periods[periodsCompleted],dt_frac_backtrack,sumPeriods);
            periodsCompleted++;

            // Find temperature components
            timeTemp[periodsCompleted] = ((float)(itGlobal) + dt_frac_backtrack)*dt;
            cudaMemcpy(vx, d_vx, (*np)*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(vr, d_vr, (*np)*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(vt, d_vt, (*np)*sizeof(float), cudaMemcpyDeviceToHost);
            FindTemperature(tempX,tempR,tempT,vx,vr,vt,periodsCompleted,(*np), fusionVelocity);

            // CREATE PARTICLES
            #if CUDATIMERS
            cudaEventRecord(startCreate);
            startCreateCPU = clock();
            #endif
            
            int npAdd = npBuffer - (*np);
            int npOld = (*np);
            (*np) += npAdd;
            for (int p = 0; p < npAdd; p++){
                float pxRand = 0;
                cudaMemcpy(d_px+npOld+p, &pxRand, sizeof(float), cudaMemcpyHostToDevice);
                float prRand = creationRadius*sqrtf((float)rand()/(float)(RAND_MAX));
                cudaMemcpy(d_pr+npOld+p, &prRand, sizeof(float), cudaMemcpyHostToDevice);
                float vxZero = fusionVelocity;
                cudaMemcpy(d_vx+npOld+p, &vxZero, sizeof(float), cudaMemcpyHostToDevice);
                float vrZero = 0.0f;
                cudaMemcpy(d_vr+npOld+p, &vrZero, sizeof(float), cudaMemcpyHostToDevice);
                float vtZero = 0.0f;
                cudaMemcpy(d_vt+npOld+p, &vtZero, sizeof(float), cudaMemcpyHostToDevice);
                int npSet = periodsCompleted;  // the newly created particles should act as though they completed this many periods so it doesn't throw off our period calculations
                cudaMemcpy(d_pP+npOld+p, &npSet, sizeof(int), cudaMemcpyHostToDevice);
            }
            
             // fix energy (conservation of energy)
            if (energyBump){
                // Find total energy of particles
                findTotalEnergy<<<blocks, threads>>>(d_E,d_px,d_pr,d_vx,d_vr,d_vt,d_phi,dx,dr,nx,nr,nxBound,nrBound,nInside,xBound,c,C,drSlope,(*np),q,m,macroWeight);
                findKineticEnergy<<<blocks, threads>>>(d_vx, d_vr, d_vt,d_KE,m,(*np),macroWeight);
                cudaMemcpy(E, d_E, (*np)*sizeof(float), cudaMemcpyDeviceToHost);
                cudaMemcpy(KE, d_KE, (*np)*sizeof(float), cudaMemcpyDeviceToHost);
                Etotal = 0.0f;
                KEtotal = 0.0f;
                for(int p = 0; p < (*np); p++) Etotal += E[p];
                for(int p = 0; p < (*np); p++) KEtotal += KE[p];
                // bump it down to the initial energy
                float bumpRatio = sqrtf(1.0f + (Einitial-Etotal)/KEtotal);
                printf("bumpRatio = %f, Einitial = %e, Etotal = %e, KEtotal = %e\n",bumpRatio,Einitial,Etotal,KEtotal);
                bumpParticleEnergy<<<blocks, threads>>>(bumpRatio, d_vx, d_vr, d_vt, (*np));
            }
            
            # if CUDATIMERS
            cudaEventRecord(stopCreate);
            cudaEventSynchronize(stopCreate);
            cudaEventElapsedTime(&tempTime, startCreate, stopCreate);
            timeCreate += tempTime;
            timeCreateCPU += ((float)(clock() - startCreateCPU)) / CLOCKS_PER_SEC*1000;
            #endif
        }
        lastMeanP = meanP;
        lastMeanPfloat = meanPfloat;

        it++;
        itGlobal++;
    }
    
    for(int i = it+1; i < nt; i++){  // set to nan
        fusionRateTotal[i] = -1.0f;
        energyLossTotal[i] = -1.0f;
    }

    findPotentialEnergy<<<blocks, threads>>>(
                d_px, d_pr,
                d_phiVec, d_PE,
                dx, dr,
                nx, nr, 
                nxBound, nrBound, nInside, xBound,
                c, C, drSlope,
                (*np), q, macroWeight);
    
    findKineticEnergy<<<blocks, threads>>>(
                d_vx, d_vr, d_vt,
                d_KE,
                m,
                (*np),macroWeight);
    
//     printf("finished loop\n");
    Timers[3] = ((float)(clock() - startLoop)) / CLOCKS_PER_SEC*1000;
    #if CUDATIMERS
    Timers[4] = timeCreate;
    Timers[5] = timeDensity;
    Timers[6] = timePotential;
    Timers[7] = timeEfield;
    Timers[8] = timeAcceleration;
    Timers[9] = timeCollisions;
    Timers[10] = timeMove;
    Timers[11] = timeBoolTransfer;
    Timers[12] = timeKill;
    #endif
//     printf("final np = %i\n",(*np));
//     
//     %  1 - Total time including async Matlab function call overhead
//     %  2 - Total time including Matlab wrapper function call overhead
//     %  3 - Total time including CUDA function call overhead
//     %  4 - Total time of mexCUDA loop
//     %  5 - InLoop: Create
//     %  6 - InLoop: Density
//     %  7 - InLoop: Potential
//     %  8 - InLoop: Electric field
//     %  9 - InLoop: Acceleration
//     % 10 - InLoop: Collisions
//     % 11 - InLoop: Move
//     % 12 - InLoop: BoolTransfer
//     % 13 - InLoop: Kill
//     % 14 - Total InLoop
//     % 15 - CUDA overhead
//     % 16 - Matlab wrapper function
//     % 17 - Matlab async overhead (+plotting if applicable)
// 
      
    // wait for the kernel to finish
    //cudaDeviceSymnchronize(); 
    // copy output vector from GPU buffer to host memory
//     printf("copying necessary data from GPU to CPU\n");
    cudaMemcpy(px, d_px, npBuffer*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(pr, d_pr, npBuffer*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vx, d_vx, npBuffer*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vr, d_vr, npBuffer*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vt, d_vt, npBuffer*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(rho, d_rho, nPM*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(phi, d_phi, nPM*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(Ex, d_Ex, nEM*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(Er, d_Er, nEM*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(phiVec, d_phiVec, nEM*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(fusionRate, d_fusionRate, nPM*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cellResidents, d_cellResidents, maxParticlesPerCell*nPM*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(partNs, d_partNs, maxParticlesPerCell*nPM*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(partAs, d_partAs, maxParticlesPerCell*nPM*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(PE, d_PE, npBuffer*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(KE, d_KE, npBuffer*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cellOccupants, d_cellOccupantsMax, nPM*sizeof(float), cudaMemcpyDeviceToHost);
    int maxOcc = 0;
    int cellAddr = 0;
    for (int cc = 0; cc < nPM; cc++){
        if (cellOccupants[cc] > maxOcc){
            maxOcc = cellOccupants[cc];
            cellAddr = cc;
        }
    }
    if(maxOcc > maxParticlesPerCell){
        printf("Warning!  you had %i particles in cell %i, and your maximum allocation is %i!!  increase maxPartilcesPerCell\n",maxOcc,cellAddr,maxParticlesPerCell);
    }else{
        printf("highest collision cell occupancy was in cell %i.... %i of %i spots occupied\n",cellAddr,maxOcc,maxParticlesPerCell);
    }
     
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    
    printf("entered mex function\n");
    clock_t startCUDA = clock();
     

    
    int i = 0;
    // Get scalars and pointers of input data
    int      load                = (int   )  mxGetScalar(prhs[i]); i++;
    int     *np                  = (int*  )  mxGetData  (prhs[i]); i++;
    float   *remainder           = (float*)  mxGetData  (prhs[i]); i++;
    float   *periods             = (float*)  mxGetData  (prhs[i]); i++;
    float   *px                  = (float*)  mxGetData  (prhs[i]); const mwSize npBuffer =  mxGetNumberOfElements (prhs[i]); i++;
    float   *pr                  = (float*)  mxGetData  (prhs[i]); i++;
    float   *vx                  = (float*)  mxGetData  (prhs[i]); i++;
    float   *vr                  = (float*)  mxGetData  (prhs[i]); i++;
    float   *vt                  = (float*)  mxGetData  (prhs[i]); i++;
    float   *KE                  = (float*)  mxGetData  (prhs[i]); i++;
    float   *PE                  = (float*)  mxGetData  (prhs[i]); i++;
    float   *Ex                  = (float*)  mxGetData  (prhs[i]); const mwSize nEM      =  mxGetNumberOfElements (prhs[i]); i++;
    float   *Er                  = (float*)  mxGetData  (prhs[i]); i++;
    float   *Bx                  = (float*)  mxGetData  (prhs[i]); i++;
    float   *Br                  = (float*)  mxGetData  (prhs[i]); i++;
    float   *rho                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *phi                 = (float*)  mxGetData  (prhs[i]); const mwSize nPM      =  mxGetNumberOfElements (prhs[i]); i++;
    float   *phiVec              = (float*)  mxGetData  (prhs[i]); i++;
    float    dx                  =           mxGetScalar(prhs[i]); i++;
    float    dr                  =           mxGetScalar(prhs[i]); i++;
    int      nx                  = (int)     mxGetScalar(prhs[i]); i++;
    int      nr                  = (int)     mxGetScalar(prhs[i]); i++;
    int      nxBound             = (int)     mxGetScalar(prhs[i]); i++;
    int      nrBound             = (int)     mxGetScalar(prhs[i]); i++;
    int      nInside             = (int)     mxGetScalar(prhs[i]); i++;
    float    xBound              = (float)   mxGetScalar(prhs[i]); i++;
    float    zOff                = (float)   mxGetScalar(prhs[i]); i++;
    float    c                   = (float)   mxGetScalar(prhs[i]); i++;
    float    C                   = (float)   mxGetScalar(prhs[i]); i++;
    float    drSlope             = (float)   mxGetScalar(prhs[i]); i++;
    float    dt                  =           mxGetScalar(prhs[i]); i++;
    mwSize   nt                  = (mwSize)  mxGetScalar(prhs[i]); i++;
    mwSize   numPeriods          = (mwSize)  mxGetScalar(prhs[i]); i++;
    float    creationRate        = (float)   mxGetScalar(prhs[i]); i++;
    float    beamPoint           = (float)   mxGetScalar(prhs[i]); i++;
    float    creationRadius      = (float)   mxGetScalar(prhs[i]); i++;
    float    fusionVelocity      = (float)   mxGetScalar(prhs[i]); i++;
    float    LxZero              = (float)   mxGetScalar(prhs[i]); i++;
    float    LrZero              = (float)   mxGetScalar(prhs[i]); i++;
    float    tanangle            = (float)   mxGetScalar(prhs[i]); i++;
    float    q                   = (float)   mxGetScalar(prhs[i]); i++;
    float    m                   = (float)   mxGetScalar(prhs[i]); i++;
    float    qom                 = (float)   mxGetScalar(prhs[i]); i++;
    float   *PM1                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *PM2                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *PM3                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *PM4                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *PM5                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *PM6                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *PM7                 = (float*)  mxGetData  (prhs[i]); i++;
    int      PM_OFF_1A           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_1B           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_1T           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_2A           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_2B           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_2T           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_6A           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_6B           = (int)     mxGetScalar(prhs[i]); i++; 
    int      PM_OFF_6T           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_7A           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_7B           = (int)     mxGetScalar(prhs[i]); i++;
    int      PM_OFF_7T           = (int)     mxGetScalar(prhs[i]); i++;
    float   *EX1                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *EX2                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *EX3                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *EX4                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *EX5                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *ER1                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *ER2                 = (float*)  mxGetData  (prhs[i]); i++;
    float   *ER3                 = (float*)  mxGetData  (prhs[i]); i++;
    int      EX_OFF_1A           = (int)     mxGetScalar(prhs[i]); i++;
    int      EX_OFF_1B           = (int)     mxGetScalar(prhs[i]); i++;
    int      EX_OFF_1T           = (int)     mxGetScalar(prhs[i]); i++;
    int      EX_OFF_5A           = (int)     mxGetScalar(prhs[i]); i++;
    int      EX_OFF_5B           = (int)     mxGetScalar(prhs[i]); i++;
    int      EX_OFF_5T           = (int)     mxGetScalar(prhs[i]); i++;
    int      ER_OFF_1A           = (int)     mxGetScalar(prhs[i]); i++;
    int      ER_OFF_1B           = (int)     mxGetScalar(prhs[i]); i++;
    int      ER_OFF_1T           = (int)     mxGetScalar(prhs[i]); i++;
    int      ER_OFF_3A           = (int)     mxGetScalar(prhs[i]); i++;
    int      ER_OFF_3B           = (int)     mxGetScalar(prhs[i]); i++;
    int      ER_OFF_3T           = (int)     mxGetScalar(prhs[i]); i++;
    float   *CellVolumesPM       = (float*)  mxGetData  (prhs[i]); i++;
    float   *bmaxPM              = (float*)  mxGetData  (prhs[i]); i++;
    float    a_coeff             = (float)   mxGetScalar(prhs[i]); i++;
    float    macroWeight         = (float)   mxGetScalar(prhs[i]); i++;
    float    W                   = (float)   mxGetScalar(prhs[i]); i++;
    int      itPM                = (int)     mxGetScalar(prhs[i]); i++;
    float    L                   = (float)   mxGetScalar(prhs[i]); i++;
    float   *bDirichlet          = (float*)  mxGetData  (prhs[i]); i++;
    float   *zOffPM              = (float*)  mxGetData  (prhs[i]); i++;
    int      maxParticlesPerCell = (int)     mxGetScalar(prhs[i]); i++;
    int     *cellResidents       = (int*)    mxGetData  (prhs[i]); i++;
    float   *partNs              = (float*)  mxGetData  (prhs[i]); i++;
    float   *partAs              = (float*)  mxGetData  (prhs[i]); i++;
    int     *cellOccupants       = (int*)    mxGetData  (prhs[i]); i++; 
    int      extSeed             = (int)     mxGetScalar(prhs[i]); i++;
    bool     FixEnergy           = (bool)    mxGetScalar(prhs[i]); i++;
    bool     loadPrevMode        = (bool)    mxGetScalar(prhs[i]); i++;
    bool     energyBump          = (bool)    mxGetScalar(prhs[i]); i++;
    float   *Timers              = (float*)  mxGetData  (prhs[i]); i++;
    float   *fusionRate          = (float*)  mxGetData  (prhs[i]); i++;


    
    int j = 0;
    plhs[j]  = mxCreateNumericArray(1, &nt, mxSINGLE_CLASS, mxREAL);
    float*  energyLossTotal;//  = (float*) plhs[j];
    energyLossTotal = (float*)mxGetData(plhs[j]);
    j++;
    plhs[j]  = mxCreateNumericArray(1, &nt, mxSINGLE_CLASS, mxREAL);
    float*  fusionRateTotal;//  = (float*) plhs[j];
    fusionRateTotal = (float*)mxGetData(plhs[j]);

 
    printf("launching kernel\n");
    kernelFunction(
            load,
            px, pr, vx, vr, vt, Ex, Er, Bx, Br, rho, phi, phiVec,
            KE, PE,
            dx, dr, nx, nr, nEM, nPM,
            nxBound, nrBound, nInside, xBound, zOff,
            c, C, drSlope,
            dt, nt, numPeriods, np,
            creationRate, remainder, beamPoint, creationRadius, fusionVelocity, LxZero, LrZero, tanangle, q, m, qom, npBuffer,
            PM1,PM2,PM3,PM4,PM5,PM6,PM7,
            PM_OFF_1A,PM_OFF_1B,PM_OFF_1T,PM_OFF_2A,PM_OFF_2B,PM_OFF_2T,PM_OFF_6A,PM_OFF_6B,PM_OFF_6T,PM_OFF_7A,PM_OFF_7B,PM_OFF_7T,
            EX1,EX2,EX3,EX4,EX5,ER1,ER2,ER3,
            EX_OFF_1A,EX_OFF_1B,EX_OFF_1T,EX_OFF_5A,EX_OFF_5B,EX_OFF_5T,ER_OFF_1A,ER_OFF_1B,ER_OFF_1T,ER_OFF_3A,ER_OFF_3B,ER_OFF_3T,
            CellVolumesPM,bmaxPM,a_coeff,macroWeight,W,itPM,1.0f,bDirichlet,zOffPM,
            maxParticlesPerCell,cellResidents,partNs,partAs,fusionRate,cellOccupants,extSeed,FixEnergy,loadPrevMode,energyBump,Timers,energyLossTotal,fusionRateTotal,
            periods);
    
    mwSize numPeriodsOutput = periodsCompleted+1;
    
    // Transfer temperature information to output variables
    j++;
    plhs[j]  = mxCreateNumericArray(1, &numPeriodsOutput, mxSINGLE_CLASS, mxREAL);
    float*  timeTempOutput;
    timeTempOutput = (float*)mxGetData(plhs[j]);
    j++;
    plhs[j]  = mxCreateNumericArray(1, &numPeriodsOutput, mxSINGLE_CLASS, mxREAL);
    float*  tempXOutput;
    tempXOutput = (float*)mxGetData(plhs[j]);
    j++;
    plhs[j]  = mxCreateNumericArray(1, &numPeriodsOutput, mxSINGLE_CLASS, mxREAL);
    float*  tempROutput;
    tempROutput = (float*)mxGetData(plhs[j]);
    j++;
    plhs[j]  = mxCreateNumericArray(1, &numPeriodsOutput, mxSINGLE_CLASS, mxREAL);
    float*  tempTOutput;
    tempTOutput = (float*)mxGetData(plhs[j]);
    j++;
     
    for(int i = 0; i <= periodsCompleted; i++){
        timeTempOutput[i] = timeTemp[i];
        tempXOutput[i] = tempX[i];
        tempROutput[i] = tempR[i];
        tempTOutput[i] = tempT[i];
    }
   

    Timers[2] = ((float)(clock() - startCUDA)) / CLOCKS_PER_SEC*1000;

}
