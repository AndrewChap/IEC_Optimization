__global__ void moveParticles(
        float *px, float *pr,
        float *vx, float *vr, float *vt,
        float *axE, float *arE,
        float *axB, float *arB,
        bool *hitWall,
        int *pP,
        const float LxZero, const float LrZero, const float tanangle,
        const float xBound,
        const int np, const float dt)
{
    int p = threadIdx.x  + blockDim.x * blockIdx.x;
    
    while(p < np){
        
        // update particle velocities from E&B accelerations using the boris method (https://www.particleincell.com/2011/vxb-rotation/)
        
        // first half of E-field acceleration
        float mx = vx[p] + 0.5f*axE[p]*dt;
        float mr = vr[p] + 0.5f*arE[p]*dt;
        float mt = vt[p];
        
        // "t-vector"
        float tx = 0.5f*axB[p]*dt;
        float tr = 0.5f*arB[p]*dt;
        float tt = 0.0f;
        
        // c = m - (m cross t)
        float cx = mx + mr*tt - mt*tr;
        float cr = mr + mt*tx - mx*tt;
        float ct = mt + mx*tr - mr*tx;
        
        // s = 2*t(1+t^2)
        float TwoOverOnePlus_tSquared = 2/(1 + tx*tx + tr*tr + tt*tt);
        float sx = tx*TwoOverOnePlus_tSquared;
        float sr = tr*TwoOverOnePlus_tSquared;
        float st = tt*TwoOverOnePlus_tSquared;
        
        // l = m + (c cross s)
        float lx = mx + cr*st - ct*sr;
        float lr = mr + ct*sx - cx*st;
        float lt = mt + cx*sr - cr*sx;
        
        // second half of E-field acceleration
        vx[p] = lx + 0.5f*axE[p]*dt;
        vr[p] = lr + 0.5f*arE[p]*dt;
        vt[p] = lt;
        
        // end of boris method
        
        // update particle velocities in 2D3V using the Birdsall method (from Birdsall's "plasma physics via computer simulation" textbook)
		float xprime = pr[p] + vr[p]*dt;
        float yprime = vt[p]*dt;
        float rprime = sqrtf(xprime*xprime + yprime*yprime);
        float sinalpha = yprime/rprime;
        float cosalpha = sqrtf(1-sinalpha*sinalpha);
        float vy1 = cosalpha*vr[p] + sinalpha*vt[p];
        float vz1 =-sinalpha*vr[p] + cosalpha*vt[p];
        vr[p] = vy1;
        vt[p] = vz1;
        // end of 2D3V Birdsall method
        
        // update particle positions from velocities.  Since our time-step is constant this is essentially the leapfrog method
        px[p] += vx[p]*dt;
        pr[p] += vr[p]*dt;
        
        // reflect particles if they cross neumann boundaries
        if (px[p]<0) {px[p] = -px[p]; vx[p] = -vx[p]; pP[p]++;}
        if (pr[p]<0) {pr[p] = -pr[p]; vr[p] = -vr[p];}
        
        // mark particles for deletion if they cross dirichlet boundaries
        hitWall[p] = 
                ((pr[p] > LrZero) && (px[p] < xBound)) // particles leave from the r-boundary of the extra core area
                ||
                ((pr[p] > px[p]*tanangle) && (px[p] >= xBound)) // particles hit the angled wall (wall between channels)
                || 
                (px[p] > LxZero); // particles hit the outer end bound of the beamline
        
        // check for NaNs and remove them from the simulation if found
        if(
                px[p]!=px[p] ||
                pr[p]!=pr[p] ||
                vx[p]!=vx[p] ||
                vr[p]!=vr[p] ||
                vt[p]!=vt[p]
                ){
            hitWall[p] = true;
        }
                
        p += blockDim.x*gridDim.x;
    }
}