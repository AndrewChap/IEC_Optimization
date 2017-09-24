// Collision constants (eq. 33)
#define KS1 1039911.280214015f
#define KS2 11.755304125515798f
#define KS3 0.003288895762902973f
#define KS4 10.412986187918094f
#define KS5 4.174107625086222f

#define KK1 67764870.00709805f
#define KK2 19.921269114379633f
#define KK3 0.003803156300400943f
#define KK4 0.488950502727711f
#define KK5 2.576315157618015f

#define KL1 192592149.4471608f
#define KL2 20.723107576469719f
#define KL3 0.0003164631432134557f
#define KL4 166.4520492662278f
#define KL5 6.193104073633183f

#define KH1 53073311.34754473f
#define KH2 22.606559682984777f
#define KH3 0.001720035240417612f
#define KH4 6.247721957791296f
#define KH5 1.618317287450503f

// fusion cross section equation constants
#define VG      6.872834791592918f
#define V400    0.914410132291908f
#define V642    1.158453345418433f
#define VS148   3.093739793140991e-1f
#define VQ235   2.413125231696255e-5f
#define VS400   8.361458900381056e-1f
#define VS5813  1.215129014697877f
#define VS1083  2.263864997278171f
#define VS2405  5.027327163854109f
#define VS3344  6.990179640718563f
#define VQ857   3.209266479482268e-2f
#define VQ234   2.392631691928657e-1f
#define VQ138   8.321513248062195e-2f
#define VQ309   4.172161344456135e-1f
#define CV0     4.118018508437671e-26f
#define CV1     2.400000000000000e-26f
#define CV2     1.105070312500000e-26f
#define AVL     1.662407892814276e-29f
#define DV0     6.898203592814370e-26f
#define DV1     6.610000000000001e-26f
#define DV2     9.711223958333338e-26f
#define DV3     8.275015476724245e-25f
#define BV0     9.155797495917257e-28f
#define AV0     2.347466090402577e-27f
#define AV1     5.179039973767553e-28f
#define AV2     1.223970646357764e-28f
#define AV3     5.188174083068728e-28f
#define SECONDSCONVERSION      1.0e-7f

__global__ void performCollisions(
        float *vx, float *vr, float *vt,
        int *Ic,
        float *rho, const float *bmax, const float *CellVolumesPM,
        const int maxParticlesPerCell, int *cellResidents, const int *cellOccupants, int *cellOccupantsMax,
        float *fusionRate,
        float *partNs, float *partAs,
        curandState *myCurandstate, 
        const float PIdt, const float a_coeff,
        const int nPM, const float L){
    
    int c = threadIdx.x  + blockDim.x * blockIdx.x;
    
    while(c < nPM){

        
        int np = cellOccupants[c];
        if(np > cellOccupantsMax[c])        // just recording the max number of occupants that have been seen in this cell
            cellOccupantsMax[c] = np;
        
        int npTrunc = np <= maxParticlesPerCell ? np : maxParticlesPerCell;   // if we have more particles than the max allowed, we need to truncate our number of particles there

        //generate random permutation of particles in the cell (Fisher Yates algorithm: https://stackoverflow.com/questions/15961119/how-to-create-a-random-permutation-of-an-array0
        for (int p = npTrunc-1; p >= 0; --p){
            //generate a random number [0, p-1]
            int r = (int)(curand(&myCurandstate[c]) % (((unsigned int)(p))+1));  // curand generates unsigned int random numbers, so we need to convert to int
            
            //swap the last element with element at random index
            int temp = cellResidents[c*maxParticlesPerCell + p];
            cellResidents[c*maxParticlesPerCell + p] = cellResidents[c*maxParticlesPerCell + r];
            cellResidents[c*maxParticlesPerCell + r] = temp;
        }

        float fusionRateHere = 0.0f;
        for (int p = 0; p < npTrunc-1; p+=2){ // right now we don't collide with an odd particle if it is odd
            
            // select particles to be collided
            int p1 = cellResidents[c*maxParticlesPerCell + p];
            int p2 = cellResidents[c*maxParticlesPerCell + p+1];
            
//             partAs[c*maxParticlesPerCell + p] = 1.0f;
            // COM velocity
            float VCx = 0.5f*(vx[p1]+vx[p2]);
            float VCr = 0.5f*(vr[p1]+vr[p2]);
            float VCt = 0.5f*(vt[p1]+vt[p2]);
            // relative velocity of p1 in the COM frame
            float VRx = vx[p1] - VCx;
            float VRr = vr[p1] - VCr;
            float VRt = vt[p1] - VCt;
            // rel velocity squared and magnitude
            float VR2 = VRx*VRx + VRr*VRr + VRt*VRt;
            float VRm = sqrtf(VR2);
            // the velocity difference
            float VDx = vx[p1] - vx[p2];
            float VDr = vr[p1] - vr[p2];
            float VDt = vt[p1] - vt[p2];
            // the velocity difference squared and magnitude
            float VD2 = VDx*VDx + VDr*VDr + VDt*VDt;
            float VDm = sqrtf(VD2);
            float VD2i = 1.0f/VD2;

            
           
            // Aside into fusion rate calculation (see thesis)
            float v = VDm*SECONDSCONVERSION;
            float vInv = 1.0f/v;
            float v2 = v*v;
            float v4 = v2*v2;
            
            float SV;
            if (v < V400){
                SV = CV0 + CV1*v2 + CV2*v4 + AVL/((v2-VS148)*(v2-VS148) + VQ235);
            }else if(v < V642){
                float peak400 = (v2 - VS400);
                float peak400squared = peak400*peak400;
                SV = DV0 + DV1*peak400 - DV2*peak400*peak400 - DV3*peak400squared*peak400squared*peak400;
            }else{
                SV = BV0 + 
                        AV0/((v2 - VS5813)*(v2 - VS5813) + VQ857) +
                        AV1/((v2 - VS1083)*(v2 - VS1083) + VQ234) +
                        AV2/((v2 - VS2405)*(v2 - VS2405) + VQ138) +
                        AV3/((v2 - VS3344)*(v2 - VS3344) + VQ309);
            }
            
            float fusionContribution = rho[c]*rho[c]/(CellVolumesPM[c]*L*L*L) * // density squared multiplied by volume (Have to divide rho by volume to get density, so in this case we just divide once by volume)
                    vInv*expf(-VG*vInv)*(1.0f/SECONDSCONVERSION)*SV/4.0f/(float(npTrunc));  // later (in MATLAB plotting) it gets multiplied by dE (8.7e6 MeV)
            
           
            if (fusionContribution != fusionContribution){ //if nan
                fusionContribution = 0.0f;
            }
            fusionRateHere += fusionContribution;
            

                    
            // back to collision modeling
            float N = (((PIdt*VDm)*bmax[c]*bmax[c])*rho[c])/(CellVolumesPM[c]*(L*L*L))*(float(np)/float(npTrunc));  // the '(float(np)/float(npTrunc))' is there to offset any 'under-colliding' that would be done if we aren't operating on all the particles
//             partNs[c*maxParticlesPerCell + p] = N;
            // inverse of the magnitude of the relative velocity
            float VRi = VRm>0 ? 1.0f/VRm : 0.0f;
            
            if (N > 0){
                float theta;
                float a = a_coeff*VD2i/bmax[c];
                float sTildeAsqoN = (a*sqrtf(N))*(KS4 - KS1*expf(-KS2*powf(N,KS3))); // eq. 31a
                float s = sTildeAsqoN*pow((1.0f + pow((TWOOVERPI*sTildeAsqoN),KS5)),(1.0f/KS5)); // eq. 32a
                
                float U = curand_uniform(&myCurandstate[c]);
                float V = curand_uniform(&myCurandstate[c]);
                
                if(s > PIOVERTWO-0.2f){ // if s is close to pi/2, we are pretty much isotropic.  this happens when two particles are traveling close to the same speed
                    theta = acosf(1.9999f*U-0.9998f);  // use 1.9999 and 0.9998 instead of 2 and 1 to avoid the remote possibility of an argument equal to -1 or 1
                }else{
                    float kTilde = 1.0f - KK1*expf(-KK2*powf(N,KK3));
                    float lTilde = KL1*expf(-KL2*powf(N,KL3));
                    float hTilde = KH1*expf(-KH2*powf(N,KH3));

                    float k = kTilde*expf(KK4*powf(s,KK5));
                    float l = lTilde*expf(-KL4*powf(s,KL5));
                    float h = hTilde*expf(-KH4*powf(s,KH5));
                    
                    // if any of k, l, or h are greater than 1 or less than 0, set 
                    k = k > 1.0f ? 1.0f : k;
                    l = l > 1.0f ? 1.0f : l;
                    h = h > 1.0f ? 1.0f : h;
                    k = k < 0.0f ? 0.0f : k;
                    l = l < 0.0f ? 0.0f : l;
                    h = h < 0.0f ? 0.0f : h;
                    
                    // see Fig. 3 for how this branching works
                    if (U >= h){ // theta_transition or theta_low
                        float varSigma = cosf(s)/(sinf(s)*sinf(s));
                        float invVarSigma = 1.0f/varSigma;
                        if (U >= l){ // most common route: theta_low
                             // (with extra "abs" to prevent negative number in logarithm, though U should never be low enough to reach that point)
                            // (u_low-1)/kappa + 1 > 0, so: u_low > 1-kappa 
                            theta =
                                    (varSigma>1.0e6f)
                                    ?
                                        sqrtf(-2*invVarSigma*logf(fabsf((U-1.0f)/k+1.0f))) // taylor expansion of acos(1-x) at x=0 is sqrt(2*x).  we need this taylor expansion when varSigma > 1e6 or so, because evaluating 1+1e-8 in single results in 1
                                        :
                                        (
                                            (varSigma>80.0f)
                                            ?
                                            (acosf(1.0f + invVarSigma*logf(fabsf((U-1.0f)/k+1.0f)))) // eq. 25
                                            :
                                            (acosf(invVarSigma*logf(fabsf(expf(-varSigma) + 2.0f*sinhf(varSigma)*((U-1.0f)/k+1.0f))))) // eq. 24
                                         )
                                    ;
                        }else{ // second-most common route: theta_mid
                            float thetaHighUhigh = 2.0f*atanf(a*sqrtf(N/h));
                            float thetaLowUlow = 
                                    (varSigma>1.0e6f)
                                    ?
                                        sqrtf(-2*invVarSigma*logf(fabsf((l-1.0f)/k+1.0f))) // taylor expansion of acos(1-x) at x=0 is sqrt(2*x).  we need this taylor expansion when varSigma > 1e6 or so, because evaluating 1+1e-8 in single results in 1
                                        :
                                        (
                                            (varSigma>80.0f)
                                            ?
                                            (acosf(1.0f + invVarSigma*logf(fabsf((l-1.0f)/k+1.0f)))) // eq. 25
                                            :
                                            (acosf(invVarSigma*logf(fabsf(expf(-varSigma) + 2.0f*sinhf(varSigma)*((l-1.0f)/k+1.0f))))) // eq. 24
                                         )
                                    ;

                            theta = thetaLowUlow*powf((U/l), (logf(thetaLowUlow/thetaHighUhigh)/logf(l/h)));
                        }
                    }else{ // least common route: theta_high
                        theta = 2.0f*atanf(a*sqrtf(N/U)); // eq. 21
                    }
                }

                // rel velocity normalized
                float VRNx = VRx*VRi;
                float VRNr = VRr*VRi;
                float VRNt = VRt*VRi;
                // collision scattering angle
                float sinTheta = sinf(theta);
                float cosTheta = cosf(theta);
                // azimuthal angle of collision
                float phi = TWOPI*V;
                float sinPhi = sinf(phi);
                float cosPhi = cosf(phi);
                // create vector "R" perpendicular to v_rel_com_norm
                float Rx, Rr, Rt;
                if (VRNr!=1.0f){
                    Rx = -VRNt;
                    Rr = 0.0f;
                    Rt = VRNx;
                }else{ // avoids the very very unlikely possibility of getting the zero vector
                    Rx = 0.0f;
                    Rr = VRNx;
                    Rt = -VRNt;
                }
                // vector "U" = VRN cross R
                float Ux = VRNr*Rt - VRNt*Rr;
                float Ur = VRNt*Rx - VRNx*Rt;
                float Ut = VRNx*Rr - VRNr*Rx;
                // vector "P" = U*sin(phi) + R*cos(phi)
                float Px = Ux*sinPhi + Rx*cosPhi;
                float Pr = Ur*sinPhi + Rr*cosPhi;
                float Pt = Ut*sinPhi + Rt*cosPhi;
                // Normalize P
                float P2 = Px*Px + Pr*Pr + Pt*Pt;
                float Pi = P2>0.0f ? rsqrtf(P2) : 0.0f;
                Px *= Pi;
                Pr *= Pi;
                Pt *= Pi;
                // vector "T" is the unnormalized rel velocity rotated towards "P" by theta.  this is the new velocity of p1 in the COM frame
                float Tx = VRx*cosTheta + (Pr*VRt - Pt*VRr)*sinTheta;
                float Tr = VRr*cosTheta + (Pt*VRx - Px*VRt)*sinTheta;
                float Tt = VRt*cosTheta + (Px*VRr - Pr*VRx)*sinTheta;
                // vector "C" is the change in velocity
                float Cx = Tx - VRx;
                float Cr = Tr - VRr;
                float Ct = Tt - VRt;
                // check for NaN's and zero them out if necessary
                if (Cx!=Cx)
                    Cx = 0.0f;
                if (Cr!=Cr)
                    Cr = 0.0f;
                if (Ct!=Ct)
                    Ct = 0.0f;
                // apply collision velocity change to particles
                vx[p1] += Cx;
                vr[p1] += Cr;
                vt[p1] += Ct;
                vx[p2] -= Cx;
                vr[p2] -= Cr;
                vt[p2] -= Ct;
            }
        }
        fusionRate[c] = fusionRateHere;
//          fusionRate[c] = 1.0f;
//          partNs[c] = 1.0f;
        c += blockDim.x*gridDim.x;
    }
}

__global__ void setupRandState(int extSeed, int nPM, curandState *state){
    int idx = threadIdx.x+blockDim.x*blockIdx.x;
    while (idx < nPM){
        curand_init(extSeed, idx, 0, &state[idx]);
        idx += blockDim.x*gridDim.x;
    }
}

