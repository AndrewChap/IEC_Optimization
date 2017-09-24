__global__ void findElectricField(
        const float *phi,
        float *Ex, float *Er,
        const float *EX1,
        const float *EX2,
        const float *EX3,
        const float *EX4,
        const float *EX5,
        const float *ER1,
        const float *ER2,
        const float *ER3,
        const int EX_OFF_1A,
        const int EX_OFF_1B,
        const int EX_OFF_1T,
        const int EX_OFF_5A,
        const int EX_OFF_5B,
        const int EX_OFF_5T,
        const int ER_OFF_1A,
        const int ER_OFF_1B,
        const int ER_OFF_1T,
        const int ER_OFF_3A,
        const int ER_OFF_3B,
        const int ER_OFF_3T,
        const int nEM,
        const float L){
    
    int c = threadIdx.x  + blockDim.x * blockIdx.x;
    
    while(c < nEM){
       
        int EX_OFF1 = (c < EX_OFF_1T)*EX_OFF_1A + (c >= EX_OFF_1T)*EX_OFF_1B;
        int EX_OFF5 = (c < EX_OFF_5T)*EX_OFF_5A + (c >= EX_OFF_5T)*EX_OFF_5B;
        int ER_OFF1 = (c < ER_OFF_1T)*ER_OFF_1A + (c >= ER_OFF_1T)*ER_OFF_1B;
        int ER_OFF3 = (c < ER_OFF_3T)*ER_OFF_3A + (c >= ER_OFF_3T)*ER_OFF_3B;
        
        Ex[c] = 1/L*(
                EX1[c]*phi[c-EX_OFF1] + 
                EX2[c]*phi[c-1] + 
                EX3[c]*phi[c] + 
                EX4[c]*phi[c+1] + 
                EX5[c]*phi[c+EX_OFF5]);
        
        Er[c] = 1/L*(
                ER1[c]*phi[c-ER_OFF1] + 
                ER2[c]*phi[c] + 
                ER3[c]*phi[c+ER_OFF3]);
        
        c += blockDim.x*gridDim.x;
        
    }
}