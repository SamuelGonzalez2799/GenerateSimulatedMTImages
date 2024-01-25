function EB1Coords=getEB1Coords(proteinPositions,angleRads,varargin)
% tubCoords = getTubulinCoords(pfLengths)
% tubCoords = getTubulinCoords(pfLengths,rotationAngleRads)
% tubCoords = getTubulinCoords(pfLengths,rotationAngleRads,_Name,_Value)
%             Name value pairs of form
%               'fluorRadius' , value of distance of fluorophore from MT
%               center
%               'rotationOption' , one of 'minusEnd','center','plusEnd'
PF_NUM = 13;
FLUOR_DIST_FROM_MT_CENTER = 15;
TUB_LENGTH = 8;

if nargin>2
    for p=1:2:numel(varargin)
        switch lower(varargin{p})
            case 'fluorradius'
                FLUOR_DIST_FROM_MT_CENTER = varargin{p+1};
            case 'rotationoption'
                rotationOption = lower(varargin{p+1});
                
        end
    end

elseif nargin>1
    %do nothing
else
   angleRads=0; 
end

if(angleRads>2*pi)
    error('angleRads must be a value on the interval 0 to 2PI');
end
% if(any(size(pfLengths(:)) ~= [PF_NUM,1]))
%     error('pfLengths must be a %dx1 or 1x%d vector',PF_NUM);
% end



%randSeed=randi(intmax());
%rng(randSeed);

%{
tic;
longest=max(pfLengths);
%taggedTubBoolArray = rand(PF_NUM,longest) <=tagFrac;
taggedTubLengthArray = rand(PF_NUM,longest) <=tagFrac;
%actualTagFrac = sum(taggedTubLengthArray(:))/numel(taggedTubLengthArray)
taggedTubLengthArray = bsxfun(@times,taggedTubLengthArray,1:longest);
taggedTubLengthArray(bsxfun(@gt,taggedTubLengthArray,pfLengths(:)))=0;
taggedTubPFArray = bsxfun(@times,taggedTubLengthArray>0,[1:PF_NUM]');
taggedTubLengthArray = taggedTubLengthArray(taggedTubLengthArray(:)>0);
taggedTubPFArray = taggedTubPFArray(taggedTubPFArray(:)>0);
% x positions are the lengths of each from the origin
% y positions are calculated as the y of a circle with the first PF located
% at a y of 0.Calc using sin(). Circle radius is FLUOR_DIST_FROM_MT_CENTER
mtPoints = [TUB_LENGTH*taggedTubLengthArray , FLUOR_DIST_FROM_MT_CENTER*sin(2*pi*(taggedTubPFArray-1)/13)];
fprintf('The bsxfun method took %f seconds \n',toc);
%}

%rng(randSeed);

%tic;
nEB1Tot=size(proteinPositions,1);
EB1Coords=zeros(nEB1Tot,2);
startIdx=1;
for iPF = 1:nEB1Tot
    EB1Coords(iPF,1)=proteinPositions(iPF,2)*TUB_LENGTH;
    if proteinPositions(iPF,1)<14;
    EB1Coords(iPF,2)=FLUOR_DIST_FROM_MT_CENTER*sin(2*pi*(proteinPositions(iPF,1)-1)/13);
    else
        EB1Coords(iPF,2)=FLUOR_DIST_FROM_MT_CENTER*sin(2*pi*(proteinPositions(iPF,1)-14)/13)%for binding on fourteen edge, say its at first edge
    end
%     nTub = pfLengths(iPF);%sum(pfLengths);
%     idxArray = 1:nTub;
%     EB1Coords(startIdx:nTub+startIdx-1,1)=idxArray*TUB_LENGTH;
%     EB1Coords(startIdx:nTub+startIdx-1,2)=FLUOR_DIST_FROM_MT_CENTER*sin(2*pi*(iPF-1)/13);
%     startIdx = nTub+startIdx;
end
%EB1Coords(startIdx:end,:)=[];

%fprintf('The loop method took %f seconds \n',toc);

%Alternate Method
%create vector with lengths and pfNums
%{
lengthsVect = zeros(sum(pfLengths),1);
pfNumVect = lengthsVect
for iPF = 1:PF_NUM
lengthsVect(x:y)=1:pflenngths(iPF);
pfVect=(x:y)=iPF
end
%Then make taggVsNoTag boolean array, and take out elements from the two
created vectors.
%}

%rotation

if(angleRads~=0)
    switch rotationOption
        case 'minusend'
            rotMatrix = [cos(angleRads),sin(angleRads); -sin(angleRads),cos(angleRads)];
            EB1Coords = EB1Coords*rotMatrix;
        case 'center'
            maxVal = max(EB1Coords(:,1));
            EB1Coords(:,1) = EB1Coords(:,1) - maxVal/2;
            rotMatrix = [cos(angleRads),sin(angleRads); -sin(angleRads),cos(angleRads)];
            EB1Coords = EB1Coords*rotMatrix;
            EB1Coords(:,1) = EB1Coords(:,1) + maxVal/2;
        case 'plusend'
            maxVal = max(EB1Coords(:,1));
            EB1Coords(:,1) = EB1Coords(:,1) - maxVal;
            rotMatrix = [cos(angleRads),sin(angleRads); -sin(angleRads),cos(angleRads)];
            EB1Coords = EB1Coords*rotMatrix;
            EB1Coords(:,1) = EB1Coords(:,1) + maxVal;
    end
    %mtPoints(:,1)=mtPoints(:,1).*cos(angleRads)-mtPoints(:,2).*sin(angleRads);
    %mtPoints(:,2)=mtPoints(:,2).*cos(angleRads)+mtPoints(:,1).*sin(angleRads);
end

