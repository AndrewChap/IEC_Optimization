function B = MagFromDipoles( ptsDip, dip, dq, ptsB, fileString)
    
    % Drew Chap 5/29/2017
    % B = MagFromDipoles( ptsDip, dip, dq, ptsB, fileString)
    % Calculates the magnetic field 'B' at the locations 'ptsB' from dipole vectors 'dip' located
    % at points 'ptsDip'.  dq is the list of the nearest neighbor distance
    % from each dipole

%     fullFileName = [mfilename,fileString,'avgd','.mat'];
        fullFileName = [mfilename,fileString,'.mat'];

    if exist(fullFileName,'file') && ~isempty(fileString)
        load(fullFileName)
        disp('Potential point spheres already exist for this configuration!  loading...')
    else

    np = size(ptsB,1);
    B = zeros(np,3);
    dq2 = dq.^2;
    str = 'Magfield point %06i of %06i';
    deleteAmount = strlength(sprintf(str,0,0));
    strDelete = repmat('\b',1,deleteAmount);
    fprintf(str,0,np);
    
    for p = 1:np
        fprintf([strDelete,str],p,np);
        R = bsxfun(@minus,ptsDip,ptsB(p,:));
        R2 = max([sum(R.^2,2),dq2],[],2);
        RmNeg5 = R2.^(-2.5);
        RmNeg3 = R2.^(-1.5);
        DdotR = sum(dip.*R,2);
        B(p,:) = sum(3*bsxfun(@times,R,DdotR.*RmNeg5) - bsxfun(@times,dip,RmNeg3),1);
    end
        save(fullFileName,'B')
    end

end