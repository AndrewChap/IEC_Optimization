function [ptsDip, dip, dipVol, streamLineStartPoints,fileString] = CreateMagneticDipoles(dpm,MagnetThickness,inner2outer,Fpent,Fhex,V,PH,PV,PF,distBetweenPentAndVert)
    
    % creates normalized IEC magnet dipoles.  output 'dip' must be
    % multiplied by the Remanence of the magnets to give the correct dipole
    % strength

     disp(['in ',mfilename])

    fileString = sprintf('d=%1.3f_',dpm);
    fileString = sprintf([fileString,'t=%1.3f_'],MagnetThickness);
    fileString = sprintf([fileString,'I=%1.3f_'],inner2outer);
    fileString = sprintf([fileString,'O=%1.3f_'],1);
    
    fullFileName = [mfilename,fileString,'.mat'];
    if exist(fullFileName,'file')
        load(fullFileName)
        disp('Dipoles already exist for this configuration!  loading...')
    else
        
        numMagPoints = 0;
        numStreamPoints = 0;
        for i = 1:12  % loop through each pentagonal face
            PentPoint = Fpent(i,:);
            for j = 1:5   % loop through the 5 hexagons that border this pentagon
                % find the two vertices that are close to the faces
                HexPoint = Fhex(PH(i,j),:);
                for k = 1:5 % loop through 5 pentagon vertices to find the two are join this pentagon with hexagon
                    if (norm(V(PV(i,k),:) - HexPoint) < distBetweenPentAndVert)
                        VertPoint = V(PV(i,k),:);
                        TopPoint = .5*(V(PV(i,k),:) + V(PF(i,k),:));
                        TopPoint(abs(TopPoint)<0.001) = 0;
                        for m = [1:k-1, k+1:5] % loop through the other 4 vertices to find it's match
                            if (norm(V(PV(i,m),:) - HexPoint) < distBetweenPentAndVert)
                                HingePoint = .5*(V(PV(i,m),:) + VertPoint);
                                HingePoint(abs(HingePoint)<0.001) = 0;
                            end
                        end
                        
                    else
                        continue
                    end
                    
                    MagPointsHere = SetMagneticPoints(dpm,...
                        VertPoint,...
                        TopPoint,...
                        HingePoint);
                    
                    numMagPointsCreated = size(MagPointsHere,1);
                    ptsDipRedundant((numMagPoints+(1:numMagPointsCreated)),:) = MagPointsHere;
                    numMagPoints = numMagPoints + numMagPointsCreated;
                    %
                    %                  ChargePointsHere = SetPotentialPoints(PentPoint,...
                    %                 HexPoint,...
                    %                 VertPoint,...
                    %                 TopPoint,...
                    %                 HingePoint,...
                    %                 WallThickness,...
                    %                 1);
                    
                    
                    StreamPointsHere = SetMagneticStreamlineStartPoints(dpm,...
                        PentPoint,...
                        HexPoint,...
                        VertPoint,...
                        TopPoint,...
                        HingePoint,...
                        1.2*MagnetThickness,...
                        InnerRadius);
                    
                    numStreamPointsCreated = size(StreamPointsHere,1);
                    StreamPoints((numStreamPoints+(1:numStreamPointsCreated)),:) = StreamPointsHere;
                    numStreamPoints = numStreamPoints + numStreamPointsCreated;
                end
            end
        end
        
        
        ptsDipNorm = unique(ptsDipRedundant,'rows');
        dp = nnsearchDrew(ptsDipNorm);
        
        dxMarg = dpm;
        
        nM = ceil(sqrt(OuterRadius)*(sqrt(OuterRadius)-sqrt(InnerRadius))/dxMarg);  % not sure exactly why this line works, but it does
        
        nMsub = nM*2+1;
        mSpace = linspace(sqrt(InnerRadius),sqrt(OuterRadius),nMsub).^2;
        magCenters = mSpace(2:2:end);
        dm = mSpace(3:2:end)-mSpace(1:2:end-2);
        
        nP = size(ptsDipNorm,1);
        ptsDip = zeros(nP*nM,3);
        dpAll = zeros(nP*nM,1);
        dmAll = zeros(nP*nM,1);
        magR = zeros(nP*nM,1);
        
        for m = 1:nM
            
            range = (1:nP) + (m-1)*nP;
            ptsDip(range,:) = ptsDipNorm*magCenters(m);
            dpAll(range) = dp.*magCenters(m);
            dmAll(range) = dm(m);
            magR(range) = magCenters(m);
            
        end
        
        %     magR2 = normDrew(ptsDip);
        %     dpAll2 = dp*magCenters;
        %     dpAll2 = dpAll2(:);
        
        magR = sum(ptsDip.^2,2).^.5; % distance of magnet center from origin
        
        
        %     keyboard
        %length = dpAll, width = magR*2*MagnetThickness
        dipVol = dpAll.*(magR*2*MagnetThickness).*dmAll;
        dip = Normalize(ptsDip).*dipVol;
        dipMag = normDrew(dip);
        %     for m = 1:nM
        %
        %         range = (1:nP) + (m-1)*nP;
        %         dipTotal(m) = sum(dipMag(range));
        %         surface(m) = 4*pi*magCenters(m)^2;
        %     end
        %
        %     figure(324)
        %     plot(magCenters,dipTotal./surface,'-o')
        %     grid on
        %     plot(magCenters,dipTotal./surface,'-o',magCenters(1,end),dipTotal(1,end)./surface(1,end),'-o')
        
        streamLineStartPoints = unique(StreamPoints,'rows');
        streamLineStartPoints(2:2:end) = streamLineStartPoints(2:2:end)*1.3;
        save(fullFileName,'ptsDip', 'dip', 'dipVol', 'streamLineStartPoints')
    end
end