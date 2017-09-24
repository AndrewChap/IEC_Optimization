
__global__ void findPotentialEnergy(
        float *px, float *pr,
        float *phi, float *PE,
        const float dx, const float dr,
        const int nx, const int nr,
        const int nxBound, const int nrBound, const int nInside, const float xBound,
        const float c, const float C, const float drSlope,
        const int np, const float q, const float macroWeight){
    
    int p = threadIdx.x  + blockDim.x * blockIdx.x;
    while(p < np){
       
        float fix;
        float drHere;
        float xShape, xLeft, xRight, ix;
        
        if (px[p] > xBound){
            fix = ((float)(nxBound-1))+(sqrtf((px[p]-xBound)*C+1.0f)-1.0f)/c;  // use nxBound instead of nxBoundPM here
            ix = floor(fix);
            xLeft  = xBound + ((c*(ix-(float)nxBound + 1.0f) + 1.0f)*(c*(ix-(float)nxBound + 1.0f) + 1.0f) - 1.0f)/C;
            xRight = xBound + ((c*(ix-(float)nxBound + 2.0f) + 1.0f)*(c*(ix-(float)nxBound + 2.0f) + 1.0f) - 1.0f)/C;
            xShape = (px[p]-xLeft)/(xRight-xLeft);
            drHere = dr + (px[p]-xBound)*drSlope;
        }else{
            fix = px[p]/dx;
            drHere = dr;
            ix = floor(fix);
            xShape = fix - ix;
        }
        
        float fir = pr[p]/drHere;
        float ir = floor(fir);
        float rminus = drHere*ir;
        float rplus = rminus + drHere;
        // Equation 4.2 in Ruyten "Density-Conserving Shape Factors for Particle Simulations in Cylindrical and Spherical Coordinates" (1991)
//         float rShape = (rplus - pr[p])*(2*rplus + 3*rminus - pr[p])/(2*rplus*rplus - rminus*rminus);
        // Equation 4.3 gives better results!!!
        float rShape = (rplus - pr[p])/(rplus-rminus)*(3.0f + rminus/pr[p])/4.0f;
        
        int IX = (int)ix;
        int IR = (int)ir;
        
        int Imm, Ipm, Imp, Ipp; // The four indices to where we are interpolating our density for this particle
        
        IndexMapper(&Imm, &Ipm, &Imp, &Ipp, IX, IR, nxBound, nrBound, nx, nInside);
        
        PE[p] =  macroWeight*q*
                    ((1.0f-xShape)*(     rShape)*phi[Imm]
                  +  (     xShape)*(     rShape)*phi[Ipm]
                  +  (1.0f-xShape)*(1.0f-rShape)*phi[Imp]
                  +  (     xShape)*(1.0f-rShape)*phi[Ipp]);
        
        p += blockDim.x*gridDim.x;
    }
}

__global__ void findKineticEnergy(
        float *vx, float *vr, float *vt,
        float *KE,
        const float m,
        const int np, const float macroWeight){
    
    int p = threadIdx.x  + blockDim.x * blockIdx.x;
    while(p < np){
        
        KE[p] = macroWeight*0.5f*m*(vx[p]*vx[p] + vr[p]*vr[p] + vt[p]*vt[p]);
        
        p += blockDim.x*gridDim.x;
    }
}   

// __global__ void sumEnergy(
//         float *KE, float *PE, float *E){
//     
//     int p = threadIdx.x  + blockDim.x * blockIdx.x;
//     while(p < np){
//         
//         E[p] = KE[p] + PE[p];
//         
//         p += blockDim.x*gridDim.x;
//     }
// }   

__global__ void findTotalEnergy(
        float *E,
        float *px, float *pr, float *vx, float *vr, float *vt,
        float *phi,
        const float dx, const float dr,
        const int nx, const int nr,
        const int nxBound, const int nrBound, const int nInside, const float xBound,
        const float c, const float C, const float drSlope,
        const int np, const float q, const float m, const float macroWeight
){

    int p = threadIdx.x  + blockDim.x * blockIdx.x;
    while(p < np){
       
        float fix;
        float drHere;
        float xShape, xLeft, xRight, ix;
        
        if (px[p] > xBound){
            fix = ((float)(nxBound-1))+(sqrtf((px[p]-xBound)*C+1.0f)-1.0f)/c;  // use nxBound instead of nxBoundPM here
            ix = floor(fix);
            xLeft  = xBound + ((c*(ix-(float)nxBound + 1.0f) + 1.0f)*(c*(ix-(float)nxBound + 1.0f) + 1.0f) - 1.0f)/C;
            xRight = xBound + ((c*(ix-(float)nxBound + 2.0f) + 1.0f)*(c*(ix-(float)nxBound + 2.0f) + 1.0f) - 1.0f)/C;
            xShape = (px[p]-xLeft)/(xRight-xLeft);
            drHere = dr + (px[p]-xBound)*drSlope;
        }else{
            fix = px[p]/dx;
            drHere = dr;
            ix = floor(fix);
            xShape = fix - ix;
        }
        if (xShape != xShape) xShape = 0.0f;  // Every once in a longtime we get a NaN here at x=0.  not sure why, so fixing it here
        float fir = pr[p]/drHere;
        float ir = floor(fir);
        float rminus = drHere*ir;
        float rplus = rminus + drHere;
        // Equation 4.2 in Ruyten "Density-Conserving Shape Factors for Particle Simulations in Cylindrical and Spherical Coordinates" (1991)
//         float rShape = (rplus - pr[p])*(2*rplus + 3*rminus - pr[p])/(2*rplus*rplus - rminus*rminus);
        // Equation 4.3 gives better results!!!
        float rShape = (rplus - pr[p])/(rplus-rminus)*(3.0f + rminus/pr[p])/4.0f;
        if (rShape != rShape) rShape = 0.0f;   // I don't think this ever happens, but since it can happen with xShape might as well put it here to be safe
        int IX = (int)ix;
        int IR = (int)ir;
        
        int Imm, Ipm, Imp, Ipp; // The four indices to where we are interpolating our density for this particle
        
        IndexMapper(&Imm, &Ipm, &Imp, &Ipp, IX, IR, nxBound, nrBound, nx, nInside);
        
        float PE =  macroWeight*q*
                    ((1.0f-xShape)*(     rShape)*phi[Imm]
                  +  (     xShape)*(     rShape)*phi[Ipm]
                  +  (1.0f-xShape)*(1.0f-rShape)*phi[Imp]
                  +  (     xShape)*(1.0f-rShape)*phi[Ipp]);
        
//         if (PE != PE) PE = 0.0f;
        
        float KE = macroWeight*0.5f*m*(vx[p]*vx[p] + vr[p]*vr[p] + vt[p]*vt[p]);
        E[p] = PE + KE;
        
        p += blockDim.x*gridDim.x;
    }
}

// change the total energy by scaling the velocity of each particle
__global__ void bumpParticleEnergy(
        const float bumpRatio,
        float *vx, float *vr, float *vt,
        const int np
){
    int p = threadIdx.x  + blockDim.x * blockIdx.x;
    while(p < np){
        vx[p] *=bumpRatio;
        vr[p] *=bumpRatio;
        vt[p] *=bumpRatio;
        p += blockDim.x*gridDim.x;
    }
}