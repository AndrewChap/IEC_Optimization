__device__ void IndexMapper(
        int *Imm,
        int *Ipm,
        int *Imp,
        int *Ipp,
        const int IX,
        const int IR,
        const int nxBound,
        const int nrBound,
        const int nx,
        const int nInside
        ){
    
    // maps our x and r indices to the linear indices of our computational mesh.  IX and IR correspond to Imm, and then Ipm, Imp, and Ipp are the 3 offsets that complete the quadrangle
    if (IR > nrBound){
        Imm[0] =  (nInside + IX + (IR - nrBound)*nxBound);
        Ipm[0] =  (Imm[0] + 1);
        Imp[0] =  (nInside + IX + (IR - nrBound+1)*nxBound);
        Ipp[0] =  (Imp[0] + 1);
    }else if(IR == nrBound){
        Imm[0] =  (IX + IR*nx);
        Ipm[0] =  (Imm[0] + 1);
        Imp[0] =  (Imm[0] + nxBound);
        Ipp[0] =  (Imp[0] + 1);
    }else{
        Imm[0] =  (IX + IR*nx);
        Ipm[0] =  (Imm[0] + 1);
        Imp[0] =  (Imm[0] + nx);
        Ipp[0] =  (Imp[0] + 1);
    }
    
}