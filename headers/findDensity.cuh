__global__ void findDensity(    // finds the number of particles in each cell
        float *px, float *pr,
        int *Ic,
        float *rho,
        const float dx, const float dr,
        const int nxPM, const int nrPM,
        const int nxBoundPM, const int nrBoundPM, const float xBound,
        const float c, const float C, const float drSlope, const float macroWeight,
        const int maxParticlesPerCell,  int *cellResidents, int *cellOccupants,
        const int np){
    
    int p = threadIdx.x  + blockDim.x * blockIdx.x;
    int nInsidePM = nxPM*nrBoundPM;
    
    while(p < np){
        
        float fix;
        float drHere;
        float xShape, xLeft, xRight, ix;
        
        if (px[p] > xBound){
            fix = ((float)(nxBoundPM))+(sqrtf((px[p]-xBound)*C+1.0f)-1.0f)/c;  // use nxBound instead of nxBoundPM here
            ix = floor(fix);
            xLeft  = xBound + ((c*(ix-(float)nxBoundPM     ) + 1)*(c*(ix-(float)nxBoundPM     ) + 1) - 1)/C;
            xRight = xBound + ((c*(ix-(float)nxBoundPM+1.0f) + 1)*(c*(ix-(float)nxBoundPM+1.0f) + 1) - 1)/C;
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
        float rShape = (rplus - pr[p])/drHere*(3.0f + rminus/pr[p])/4.0f;
        if (rShape != rShape) rShape = 0.0f;   // I don't think this ever happens, but since it can happen with xShape might as well put it here to be safe

        int IX = (int)ix;
        int IR = (int)ir;
        
        register int Imm, Ipm, Imp, Ipp; // The four indices to where we are interpolating our density for this particle
        int cellOccOld;

        if ((IX < nxPM-1 && IR < nrPM-1) && !(IX >= nxBoundPM-1 && IR >= nrBoundPM-1 )){
        
            IndexMapper(&Imm, &Ipm, &Imp, &Ipp, IX, IR, nxBoundPM, nrBoundPM, nxPM, nInsidePM);

            // find the particle's "home" where it can collide with other particles
            int IcHere =
                    Imm*(xShape <= 0.5f)*(rShape >= 0.5f) +
                    Ipm*(xShape >  0.5f)*(rShape >= 0.5f) +
                    Imp*(xShape <= 0.5f)*(rShape <  0.5f) +
                    Ipp*(xShape >  0.5f)*(rShape <  0.5f);

            cellOccOld = atomicAdd(cellOccupants+IcHere, 1);
            if (cellOccOld < maxParticlesPerCell){
                cellResidents[IcHere*maxParticlesPerCell + cellOccOld] = p;
            }

            float Wmm = macroWeight*(1.0f-xShape)*(     rShape);
            float Wpm = macroWeight*(     xShape)*(     rShape);
            float Wmp = macroWeight*(1.0f-xShape)*(1.0f-rShape);
            float Wpp = macroWeight*(     xShape)*(1.0f-rShape);
            
            if (Wmm != Wmm)
                Wmm = 0.0f;
            if (Wpm != Wpm)
                Wpm = 0.0f;
            if (Wmp != Wmp)
                Wmp = 0.0f;
            if (Wpp != Wpp)
                Wpp = 0.0f;
            
            // CUDA version of density writing
            atomicAdd(rho+Imm, Wmm);
            atomicAdd(rho+Ipm, Wpm);
            atomicAdd(rho+Imp, Wmp);
            atomicAdd(rho+Ipp, Wpp);
            
            // Serial version of density writing
            // rho[Imm] = rho[Imm] + macroWeight*(1.0f-xShape)*(     rShape));
            // rho[Ipm] = rho[Ipm] + macroWeight*(     xShape)*(     rShape));
            // rho[Imp] = rho[Imp] + macroWeight*(1.0f-xShape)*(1.0f-rShape));
            // rho[Ipp] = rho[Ipp] + macroWeight*(     xShape)*(1.0f-rShape));
            
//             // CUDA version of density writing
//             atomicAdd(rho+Imm, macroWeight*(1.0f-xShape)*(     rShape));
//             atomicAdd(rho+Ipm, macroWeight*(     xShape)*(     rShape));
//             atomicAdd(rho+Imp, macroWeight*(1.0f-xShape)*(1.0f-rShape));
//             atomicAdd(rho+Ipp, macroWeight*(     xShape)*(1.0f-rShape));
            
            
            
        }
        
        p += blockDim.x*gridDim.x;
        
    }
}