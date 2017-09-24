__global__ void findPotential(
        float *rho,
        float *phi,
        float *phiUpdate,
        float *CellVolumesPM,
        float *bDirichlet,
        float *zOffPM,
        const float *PM1,
        const float *PM2,
        const float *PM3,
        const float *PM4,
        const float *PM5,
        const float *PM6,
        const float *PM7,
        const int PM_OFF_1A,
        const int PM_OFF_1B,
        const int PM_OFF_1T,
        const int PM_OFF_2A,
        const int PM_OFF_2B,
        const int PM_OFF_2T,
        const int PM_OFF_6A,
        const int PM_OFF_6B,
        const int PM_OFF_6T,
        const int PM_OFF_7A,
        const int PM_OFF_7B,
        const int PM_OFF_7T,
        const int nPM,
        const float W,
        const float L,
        const float q){
    
    int c = threadIdx.x  + blockDim.x * blockIdx.x;
    
    while(c < nPM){
        
        int PM_OFF1 = (c < PM_OFF_1T)*PM_OFF_1A + (c >= PM_OFF_1T)*PM_OFF_1B;
        int PM_OFF2 = (c < PM_OFF_2T)*PM_OFF_2A + (c >= PM_OFF_2T)*PM_OFF_2B;
        int PM_OFF6 = (c < PM_OFF_6T)*PM_OFF_6A + (c >= PM_OFF_6T)*PM_OFF_6B;
        int PM_OFF7 = (c < PM_OFF_7T)*PM_OFF_7A + (c >= PM_OFF_7T)*PM_OFF_7B;
        
//         if(rho[c]!=rho[c]){
//             rho[c]=0.0f;  // Don't know why, but sometimes get a NaNs at the origin in rho
//         }

        phiUpdate[c] = (1-W)* phi[c] +
                          W * (PM4[c]*((PM1[c]*phi[c-PM_OFF1] + 
                                        PM2[c]*phi[c-PM_OFF2] + 
                                        PM3[c]*phi[c-1] +
                                        PM5[c]*phi[c+1] +
                                        PM6[c]*phi[c+PM_OFF6] +
                                        PM7[c]*phi[c+PM_OFF7])/(L*L) -  // poisson's equation has 1/x^2 in it
                                        q*rho[c]*zOffPM[c]/(CellVolumesPM[c]*(L*L*L))/EPS0 + bDirichlet[c]));  // Multiply volumes by L^3 to get units
        
        c += blockDim.x*gridDim.x;
        
    }
}