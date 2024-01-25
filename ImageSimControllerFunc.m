% ImageSimController
% Creates fluorophore points to pass to ImageSim
function[outputImage]=ImageSimControllerFunc(proteinPositions,protoLengths)
%% Initializations
% clear all;
% close all;

addpath('Z:Groups\LAB-klei0091\Sam\Taylor code for MT dyanmics with tpx2/MicroTubule');
PrepareFigureSize;
ImSimFunc=@ImageSimV6;%was 6 i changed it to version 8

%imageSize in x,y
imSize = [1024,32];%originally 254
%numerical aperature 100x lense
na = 1.49;
%pixel size in nm
pixSize=64;%Was 160 originally, changed to match 

%wavelength of the three color channels
lambdaR = 675;%nm (RhodamineB emission maximum is 564)
lambdaG = 510;%nm (GFP emission maximum is 509)
lambdaB = 475;%nm

%Noise and standard deviation defaults
camNoise = 1;
camStdv = 0.1;
backNoise = 10;
backStdv = 1;

%points to use for fluorophore locations
nBackFluorsPerUnitArea=1;%Originally, was 3 for 160 nm pixel size. Used 5 or 6 for 160 nm pixel size. ALso used 1 for 64 nm pixel sizze. 


%illumination mask
illuminationMask= ones(imSize(2),imSize(1));% gaussIt2D(1,imSize(1),1,imSize(2),150,imSize(1)/2,imSize(2)/2);%ones(imSize(2),imSize(1));%gaussIt2D(1,imSize(1),1,imSize(2),50,imSize(1)/2,imSize(2)/2);



%% Initialization Presets
%Presets for testing certain aspects, overwrite the above.
%preset = 'normalizationTestPoint';
preset = 'otherwise';

switch preset
    case 'rotationTest'
        imSize = [65,65];
        nBackFluorsPerUnitArea = 0;
        camNoise = 0;
        camStdv = 0;
        backNoise = 0;
        backStdv = 0;
        illuminationMask= ones(imSize(2),imSize(1));
    case 'normalizationTest'
        imSize = [64,48];
        nBackFluorsPerUnitArea = 0;
        numI = 10;
        camNoise = 0;
        camStdv = 0;
        backNoise = 0;
        backStdv = 0;
        illuminationMask= ones(imSize(2),imSize(1));
    case 'normalizationTestPoint'
        imSize = [31,31];
        nBackFluorsPerUnitArea = 0;
        numI = 16;
        camNoise = 0;
        camStdv = 0;
        backNoise = 0;
        backStdv = 0;
        illuminationMask= ones(imSize(2),imSize(1));
    otherwise
        %do nothing
end



%% Create random background fluorophores
nBackFluors = round(imSize(1)*imSize(2)*nBackFluorsPerUnitArea(1));
if(nBackFluors >0)
    pointsR = [rand(nBackFluors,1)*(imSize(1)-1)+1,rand(nBackFluors,1)*(imSize(2)-1)+1];
    pointsG=[rand(nBackFluors,1)*(imSize(1)-1)+1,rand(nBackFluors,1)*(imSize(2)-1)+1];
else
    pointsR = [];
    pointsG = [];
end

pointsB = [];


tic;
clear tubPoints
clear fluorPoints
%% Preset code
switch preset
    case 'rotationTest'
        %% Rotation Test Code
        longest=600;
        s=8;
        pfLengths=(longest-s*12:s:longest);
        rotationRads = 0;
        tubPoints = getTubulinCoords(pfLengths,rotationRads,'rotationOption','center');
        schemeValues = {'none',50,[];'poisson',50,1};
        [~,fluorIdxs] = getFluorCoords(numel(tubPoints)/2,'mixed',schemeValues);
        fluorPoints = tubPoints(fluorIdxs,:);
        fluorPoints = nm2pixels(fluorPoints,pixSize);
        numI = 5;
        %xShiftFunc = @(x) (x)*ceil((imSize(1)/(numI+1))/2)*2;
        xShiftFunc = @(x) imSize(1)/2;
        for ii=1:numI
            if ii ==1
                rotationRads(ii) = 0;
            else
                rotationRads(ii) = rand()*2*pi;
            end
            %rotate, but save the old fluorPoints
            fluorPointsNew = fluorPoints;
            rotMatrix = [cos(rotationRads(ii)),sin(rotationRads(ii)); -sin(rotationRads(ii)),cos(rotationRads(ii))];
            maxVal = max(fluorPointsNew(:,1));
            minVal = min(fluorPointsNew(:,1));
            endpoints{ii} = [minVal,min(fluorPointsNew(:,2));maxVal,max(fluorPointsNew(:,2))];%track endpoints
            fluorPointsNew(:,1) = fluorPointsNew(:,1) - maxVal/2 - minVal;
            endpoints{ii}(:,1) = endpoints{ii}(:,1) - maxVal/2 - minVal;
            fluorPointsNew = fluorPointsNew*rotMatrix;
            endpoints{ii} = endpoints{ii}*rotMatrix;
            
            %shift MT to desired location
            numMTPoints = size(fluorPointsNew,1);
            
            xShift= xShiftFunc(ii);%rand()*imSize(1);
            yShift= imSize(2)/2;%rand()*imSize(2);
            xShiftSaves(ii)=xShift;
            yShiftSaves(ii)=yShift;
            
            
            fluorPointsNew = fluorPointsNew + [zeros(numMTPoints,1)+xShift,zeros(numMTPoints,1)+yShift];
            endpoints{ii} = endpoints{ii} + [zeros(2,1)+xShift,zeros(2,1)+yShift];
            
            params = {'imSize',uint64(imSize),'lambdaR',lambdaR,'lambdaG',lambdaG,'lambdaB',lambdaB,'na',na,'pointsR',fluorPointsNew,'pointsG',pointsG,'pointsB',pointsB,'pixSize',pixSize,'camNoise',camNoise,'camStdv',camStdv,'backNoise',backNoise,'backStdv',backStdv,'illuminationMask',illuminationMask} ;
            cropMT{ii} = ImSimFunc(params);
        end
        
        buffer = xShiftFunc(1)/2;
        for kk = 1:numI
            %cropMT{kk} = outputImage(:,xShiftFunc(kk)-buffer:xShiftFunc(kk+1)-buffer, 1);
            %cropMT{kk} = outputImage(yShiftSaves(kk)+(-5:5),xShiftSaves(kk)-2 : xShiftSaves(kk)+2+max(pfLengths{kk})*8/pixSize , 1);
            cropRotateNearest{kk} = imrotate(cropMT{kk},rad2deg(rotationRads(kk)),'nearest','crop');
            cropRotateBilin{kk} = imrotate(cropMT{kk},rad2deg(rotationRads(kk)),'bilinear','crop');
            cropRotateBicub{kk} = imrotate(cropMT{kk},rad2deg(rotationRads(kk)),'bicubic','crop');
        end
        cropMTTotal = cell2mat(cropMT(:)');
        
        %analysis section
        for ii = 2:numI
            score.nearest(ii) = sum(abs(cropMT{1}(:)-cropRotateNearest{ii}(:)));
            score.bilin(ii) = sum(abs(cropMT{1}(:)-cropRotateBilin{ii}(:)));
            score.bicub(ii) = sum(abs(cropMT{1}(:)-cropRotateBicub{ii}(:)));
        end
        score.nearest = sum(score.nearest);
        score.bilin = sum(score.bilin);
        score.bicub = sum(score.bicub);
        
        figure;
        subPlot(4,1,1);
        imshow(mat2gray(cropMTTotal));
        subPlot(4,1,2);
        imshow(mat2gray(cell2mat(cropRotateNearest(:)')));
        subPlot(4,1,3);
        imshow(mat2gray(cell2mat(cropRotateBilin(:)')));
        subPlot(4,1,4);
        imshow(mat2gray(cell2mat(cropRotateBicub(:)')));
    
    case 'normalizationTest'
        %% Normalization Test Code
        longest=600;
        s=8;
        pfLengths=(longest-s*12:s:longest);
        rotationRads = 0;
        tubPoints = getTubulinCoords(pfLengths,rotationRads,'rotationOption','center');
        schemeValues = {'none',50,[];'poisson',50,1};
        [~,fluorIdxs] = getFluorCoords(numel(tubPoints)/2,'mixed',schemeValues);
        fluorPoints = tubPoints(fluorIdxs,:);
        fluorPoints = nm2pixels(fluorPoints,pixSize);
        %xShiftFunc = @(x) (x)*ceil((imSize(1)/(numI+1))/2)*2;
        xShiftFunc = @(x) imSize(1)/2;
        for ii=1:numI
            rotationRads(ii) = 0;
            %rotate, but save the old fluorPoints
            fluorPointsNew = fluorPoints;
            maxVal = max(fluorPointsNew(:,1));
            minVal = min(fluorPointsNew(:,1));
            endpoints{ii} = [minVal,min(fluorPointsNew(:,2));maxVal,max(fluorPointsNew(:,2))];%track endpoints
            fluorPointsNew(:,1) = fluorPointsNew(:,1) - maxVal/2 - minVal;
            endpoints{ii}(:,1) = endpoints{ii}(:,1) - maxVal/2 - minVal;
            
            %shift MT to desired location
            numMTPoints = size(fluorPointsNew,1);
            
            xShift= xShiftFunc(ii);%rand()*imSize(1);
            yShift= imSize(2)/2;%rand()*imSize(2);
            xShiftSaves(ii)=xShift;
            yShiftSaves(ii)=yShift;
            fluorPointsNew = fluorPointsNew + [zeros(numMTPoints,1)+xShift,zeros(numMTPoints,1)+yShift];
            endpoints{ii} = endpoints{ii} + [zeros(2,1)+xShift,zeros(2,1)+yShift];
            
            % Add increasing background sources
            nBackFluorsPerUnitArea = ii-1;
            nBackFluors = round(imSize(1)*imSize(2)*nBackFluorsPerUnitArea(1));
            if(nBackFluors >0)
                fluorPointsNew = vertcat(fluorPointsNew,[rand(nBackFluors,1)*(imSize(1)-1)+1,rand(nBackFluors,1)*(imSize(2)-1)+1]);
            end
            
            params = {'imSize',uint64(imSize),'lambdaR',lambdaR,'lambdaG',lambdaG,'lambdaB',lambdaB,'na',na,'pointsR',fluorPointsNew,'pointsG',pointsG,'pointsB',pointsB,'pixSize',pixSize,'camNoise',camNoise,'camStdv',camStdv,'backNoise',backNoise,'backStdv',backStdv,'illuminationMask',illuminationMask} ;
            cropMT{ii} = ImSimFunc(params);
            cropMT{ii} = cropMT{ii}(:,:,1);
        end
        
        bkgdRegion = {[32:48],3:size(cropMT{1},2)-2};
        [gridX,gridY] = meshgrid(bkgdRegion{2},bkgdRegion{1});
        bkgdSubs = sub2ind(size(cropMT{1}),gridY(:),gridX(:));
        cropSubBack{1} = cropMT{1};
        cropDivBack{1} = cropMT{1};
        cropSubDivBack{1} = cropMT{1};
        for kk = 2:numI
            %cropMT{kk} = outputImage(:,xShiftFunc(kk)-buffer:xShiftFunc(kk+1)-buffer, 1);
            %cropMT{kk} = outputImage(yShiftSaves(kk)+(-5:5),xShiftSaves(kk)-2 : xShiftSaves(kk)+2+max(pfLengths{kk})*8/pixSize , 1);
            cropSubBack{kk} = cropMT{kk} - median(cropMT{kk}(bkgdSubs));
            cropDivBack{kk} = cropMT{kk} / median(cropMT{kk}(bkgdSubs));
            cropSubDivBack{kk} = cropMT{kk} - median(cropMT{kk}(bkgdSubs));
            cropSubDivBack{kk} = cropSubDivBack{kk} / mean(abs(cropSubDivBack{kk}(bkgdSubs)));
        end
        cropMTTotal = cell2mat(cropMT(:)');
        
        %analysis section
        for ii = 2:numI
            score.subBack(ii) = sum(abs(cropSubBack{1}(:)-cropSubBack{ii}(:)));
            score.divBack(ii) = sum(abs(cropDivBack{1}(:)-cropDivBack{ii}(:)));
            score.subDivBack(ii) = sum(abs(cropSubDivBack{1}(:)-cropSubDivBack{ii}(:)));
        end
        score.subBack = sum(score.subBack);
        score.divBack = sum(score.divBack);
        score.subDivBack = sum(score.subDivBack);
        
        figure;
        subPlot(4,1,1);
        imshow(mat2gray(cropMTTotal));
        subPlot(4,1,2);
        imshow(mat2gray(cell2mat(cropSubBack(:)')));
        subPlot(4,1,3);
        imshow(mat2gray(cell2mat(cropDivBack(:)')));
        subPlot(4,1,4);
        imshow(mat2gray(cell2mat(cropSubDivBack(:)')));
        
    case 'normalizationTestPoint'
        %% Normalization Test Code
        fluorPoints = ones(100,2)*16;
        for ii=1:sqrt(numI)
            
            for jj = 1:sqrt(numI);
                %rotate, but save the old fluorPoints
                fluorPointsNew = fluorPoints;
                
                % Add increasing background sources
                nBackFluorsPerUnitArea = (1)/4;
                nBackFluors = round(imSize(1)*imSize(2)*nBackFluorsPerUnitArea(1));
                if(nBackFluors >0)
                    fluorPointsNew = vertcat(fluorPointsNew,[rand(nBackFluors,1)*(imSize(1)-1)+1,rand(nBackFluors,1)*(imSize(2)-1)+1]);
                end
                
                %add increasing fluorescent illumination
                illuminationFactor = jj;
                illuminationMaskNew = illuminationMask*illuminationFactor;
                
                params = {'imSize',uint64(imSize),'lambdaR',lambdaR,'lambdaG',lambdaG,'lambdaB',lambdaB,'na',na,'pointsR',fluorPointsNew,'pointsG',pointsG,'pointsB',pointsB,'pixSize',pixSize,'camNoise',camNoise,'camStdv',camStdv,'backNoise',backNoise,'backStdv',backStdv,'illuminationMask',illuminationMaskNew} ;
                cropMT{ii,jj} = ImSimFunc(params);
                cropMT{ii,jj} = cropMT{ii,jj}(:,:,1);
                bkgdRegion = {[25:30],3:size(cropMT{1},2)-2};
                [gridX,gridY] = meshgrid(bkgdRegion{2},bkgdRegion{1});
                bkgdSubs = sub2ind(size(cropMT{1}),gridY(:),gridX(:));
                allBackMed(ii,jj) = median(cropMT{ii,jj}(bkgdSubs));
            end
        end
        
        bkgdRegion = {[25:30],3:size(cropMT{1},2)-2};
        [gridX,gridY] = meshgrid(bkgdRegion{2},bkgdRegion{1});
        bkgdSubs = sub2ind(size(cropMT{1}),gridY(:),gridX(:));
        sigRegion = {[15:17],[15:17]};
        [gridX,gridY] = meshgrid(sigRegion{2},sigRegion{1});
        sigSubs = sub2ind(size(cropMT{1}),gridY(:),gridX(:));
        cropSubBack{1} = cropMT{1};
        cropDivBack{1} = cropMT{1};
        cropSubDivBack{1} = cropMT{1};
        cropDivSubBack{1} = cropMT{1};
        
        
        avgBack = mean(allBackMed(:));
        maxBack = max(allBackMed(:));
        minBack = min(allBackMed(:));
        for kk = 1:numI
            [r,c] = ind2sub([ii,jj],kk);
            %cropMT{kk} = outputImage(:,xShiftFunc(kk)-buffer:xShiftFunc(kk+1)-buffer, 1);
            %cropMT{kk} = outputImage(yShiftSaves(kk)+(-5:5),xShiftSaves(kk)-2 : xShiftSaves(kk)+2+max(pfLengths{kk})*8/pixSize , 1);
            backMed = median(cropMT{kk}(bkgdSubs));
            
            cropSubBack{kk} = cropMT{kk} - backMed;
            
            cropDivBack{kk} = cropMT{kk} / backMed;
            
            cropSubDivBack{kk} = cropMT{kk} - backMed;
            cropSubDivBack{kk} = cropSubDivBack{kk} / (backMed/mean(allBackMed(r,:)));
            
            cropDivSubBack{kk} = cropMT{kk} - backMed;
            cropDivSubBack{kk} = cropDivSubBack{kk} / backMed;
            
        end
        cropSubBack = reshape(cropSubBack,[sqrt(numI),sqrt(numI)]);
        cropDivBack = reshape(cropDivBack,[sqrt(numI),sqrt(numI)]);
        cropSubDivBack = reshape(cropSubDivBack,[sqrt(numI),sqrt(numI)]);
        cropDivSubBack = reshape(cropDivSubBack,[sqrt(numI),sqrt(numI)]);
        cropMTTotal = cell2mat(cropMT);
        maxSub = max(max(cell2mat(cropSubBack)));
        minSub = min(min(cell2mat(cropSubBack)));
        maxDiv = max(max(cell2mat(cropDivBack)));
        minDiv = min(min(cell2mat(cropDivBack)));
        maxSD = max(max(cell2mat(cropSubDivBack)));
        minSD = min(min(cell2mat(cropSubDivBack)));
        maxDS = max(max(cell2mat(cropDivSubBack)));
        minDS = min(min(cell2mat(cropDivSubBack)));
        for kk = 1:numI
            cropSubBack{kk} = mat2gray(cropSubBack{kk},[minSub,maxSub]);
            
            cropDivBack{kk} = mat2gray(cropDivBack{kk},[minDiv,maxDiv]);
            
            cropSubDivBack{kk} = mat2gray(cropSubDivBack{kk},[minSD,maxSD]);
            
            cropDivSubBack{kk} = mat2gray(cropDivSubBack{kk},[minDS,maxDS]);
        end
        
        %analysis section
        for ii = 1:numI
            score.subBack(ii) = sum(sum(abs(cropSubBack{1}(14:18,14:18)-cropSubBack{ii}(14:18,14:18))));
            score.divBack(ii) = sum(sum(abs(cropDivBack{1}(14:18,14:18)-cropDivBack{ii}(14:18,14:18))));
            score.subDivBack(ii) = sum(sum(abs(cropSubDivBack{1}(14:18,14:18)-cropSubDivBack{ii}(14:18,14:18))));
            score.divSubBack(ii) = sum(sum(abs(cropDivSubBack{1}(14:18,14:18)-cropDivSubBack{ii}(14:18,14:18))));
        end
        score.subBack = sum(score.subBack);
        score.divBack = sum(score.divBack);
        score.subDivBack = sum(score.subDivBack);
        score.divSubBack = sum(score.divSubBack);
        
        spacer = 0.005;
        figure;
        subPlot(1,1,1);
        imshow(mat2gray(cropMTTotal));
        figure;
        subPlot(1,4,1,spacer);
        imshow(mat2gray(cell2mat(cropSubBack)));
        subPlot(1,4,2,spacer);
        imshow(mat2gray(cell2mat(cropDivBack)));
        subPlot(1,4,3,spacer);
        imshow(mat2gray(cell2mat(cropSubDivBack)));
        subPlot(1,4,4,spacer);
        imshow(mat2gray(cell2mat(cropDivSubBack)));
    otherwise 
        %% Default code
        numI = 1;%adjusted by sam, was 5
        numJ = 1;
        for ii=1:numI
            for jj=1:numJ
                
                %shortest= 500;%randi(1000)+100;
                %s= ii* 8;%randi(60);
                %pfLengths=shortest:s:shortest+s*12;
                %
                longest=600;
                s=ii* 8;
                 %pfLengths(ii,jj)={(longest-s*12:s:longest)};%what it was with Taylors code
                pfLengths(ii,jj)={protoLengths};%{(longest-s*12:s:longest)};
                %
                
                % getMTFluorLocs uses protofilamentLengths, tagPercentage,rotationInRads
                % getMTFluorLocs returns locations based on the minus end being at [0,0]
                rotationRads = 0; %rand()*2*pi;
                tubPoints = getTubulinCoords(pfLengths{ii,jj},rotationRads,'rotationOption','center');
                schemeValues = {'none',85,[];'poisson',15,1};%%%%%changing poison values here....was originally 50-50
                [~,fluorIdxs] = getFluorCoords(numel(tubPoints)/2,'mixed',schemeValues);
                %fluorPoints = getMicroTubuleFluorsLocs(pfLengths,0.15,rotationRads);
                fluorPoints = tubPoints(fluorIdxs,:);
                fluorPoints = nm2pixels(fluorPoints,pixSize);
                greenPointsDimers=proteinPositions;
                greenNmPoints=getEB1Coords(greenPointsDimers,rotationRads,'rotationOption','center');%convert coordinates to distance along MT and placement along 3d circle
                greenFluorPoints=nm2pixels(greenNmPoints,pixSize);
                %shift MT to desired location
                numMTPoints = size(fluorPoints,1);
                greenNumMtPoints=size(greenFluorPoints,1);
                
                xShift= (ii-0.5) * imSize(1) / (numI+1);%rand()*imSize(1);
                yShift= imSize(2)/2;%rand()*imSize(2);
                xShiftSaves(ii,jj)=xShift;
                yShiftSaves(ii,jj)=yShift;
                
                
                fluorPoints = fluorPoints + [zeros(numMTPoints,1)+xShift,zeros(numMTPoints,1)+yShift];
                greenFluorPoints=greenFluorPoints+[zeros(greenNumMtPoints,1)+xShift,zeros(greenNumMtPoints,1)+yShift];
                if(any(fluorPoints(:)<=0))
                    %fluorophores outside bounds of image, may cause issues.
                    warning('warning');
                end
                %display(100*numel(mtPoints(:,1))/sum(pfLengths));
                pointsR = vertcat(pointsR,fluorPoints);
                pointsG=vertcat(pointsG,greenFluorPoints);
            end
        end
        
        params = {'imSize',uint64(imSize),'lambdaR',lambdaR,'lambdaG',lambdaG,'lambdaB',lambdaB,'na',na,'pointsR',pointsR,'pointsG',pointsG,'pointsB',pointsB,'pixSize',pixSize,'camNoise',camNoise,'camStdv',camStdv,'backNoise',backNoise,'backStdv',backStdv,'illuminationMask',illuminationMask} ;
        
        tic;
        outputImage = ImSimFunc(params);
        toc;
%         figure;
% %         subPlot(1,1,1);
%         imshow(mat2gray(outputImage));
end
toc;

%% End

