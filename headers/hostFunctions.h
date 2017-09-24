void FindTemperature(
        float *tempX, float *tempR, float *tempT,
        float *vx,    float *vr,    float *vt,
        const int index,
        const int np,
        const float fusionVelocity){
    // finds temperature normalized by the fusion COM energy
    
    // calculate mean velocity to subtract out of x-component
    float meanVx = 0.0f;
    for (int p = 0; p < np; p++){
        meanVx += fabs(vx[p]);
    }
    meanVx /= float(np);
    
    // calculate mean velocity to subtract out of x-component
    float meanVt = 0.0f;
    for (int p = 0; p < np; p++){
        meanVt += fabs(vt[p]);
    }
    meanVt /= float(np);

    // initialize avg square of velocity to zero
    tempX[index] = 0.0f;
    tempR[index] = 0.0f;
    tempT[index] = 0.0f;

    // caclulate avg square of velocity for each component
    for (int p = 0; p < np; p++){
        tempX[index] += (fabs(vx[p])-meanVx)*(fabs(vx[p])-meanVx);
        tempR[index] += vr[p]*vr[p];
        tempT[index] += (fabs(vt[p])-meanVt)*(fabs(vt[p])-meanVt);
    }

    // divide by number of particles to get average, and normalize by square of fusion energy
    tempX[index] /= (float(np)*fusionVelocity*fusionVelocity);
    tempR[index] /= (float(np)*fusionVelocity*fusionVelocity);
    tempT[index] /= (float(np)*fusionVelocity*fusionVelocity);

}