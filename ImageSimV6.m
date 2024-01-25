function [outputImage]=ImageSimV6(varargin)
%Simulation of fluorescence image
%   [outputI]=protFunc(parameters)
%       returns a simulated fluroscence image 'outputI'
%
%       'parameters' can be a comma seperated list of name value pairs
%        or a cell array of name value pairs. Names are NOT case sensitive
%           example:
%           protFunc('lambdaR',564,'lambdaG',509,'imSize',[512,512],...)
%           or
%           params={'lambdaR',564,'lambdaG',509,'imSize',[512,512],...}
%           protFunc(params)
%
%   Parameter list:
%   Format: ParamName(default): Description
%       lambdaR(675): emission peak of red channel (nm)
%       lambdaG(510): emission peak of green channel (nm)
%       lambdaB(475): emission peak of blue channel (nm)
%       imSize([512,512]): two element vector describing number of pixels
%       of the output image as [nPixelsInX , nPixelsInY]
%       pointsR([]): This describes the locations of red fluorophores
%           in an Nx2 array of x-y coordinates where the first column
%           is a vector of x coordinates X1...XN and the second column is
%           a vector of y coordinates Y1...YN. Positions should be in
%           reference to the desired output image and in units of pixels.
%       pointsG([]):This describes the locations of green fluorophores
%           in an Nx2 array of x-y coordinates where the first column
%           is a vector of x coordinates X1...XN and the second column is
%           a vector of y coordinates Y1...YN. Positions should be in
%           reference to the desired output image and in units of pixels.
%       pointsB([]):This describes the locations of blue fluorophores
%           in an Nx2 array of x-y coordinates where the first column
%           is a vector of x coordinates X1...XN and the second column is
%           a vector of y coordinates Y1...YN. Positions should be in
%           reference to the desired output image and in units of pixels.
%       nA(1.49): numerical aperature of theoretical microscope
%       pixSize(160): pixel size in nanometers
%       camNoise(1): mean noise value due to the camera. Units of single
%           fluorophore intensity.
%       camStdv(0.1): standard deviation of camera noise. Units of single
%           fluorophore intensity.
%       backNoise(1): mean noise value due to the background signal. Units 
%           of single fluorophore intensity. An alternative noise method is
%           to set backNoise and backStdv to 0 and include randomly
%           distributed fluorophores in your pointsR/G/B array.
%       backStdv(0.1): standard deviation of the background signal. Units 
%           of single fluorophore intensity.
%       simWindowSizeRatio(5): size of simulated window around each
%           flurophore in which it can contribute. Units are in standard
%           deviation of the fluorophore, as deterimined by its wavelength
%           and the NA of the theoretical optics, and the pixel size in
%           nanometers. Increasing simWindowSizeRatio increases accuracy 
%           with diminishing returns as you increase, at the cost of
%           increasing computation time.
%       illuminationMask(): a array the size of the image whose elements
%           indicate the relative intensity of illumination at that
%           position.Default is all ones.
%
%   Copyright Taylor Reid, Gardner Lab, University of Minnesota, 2015.
%   Please contact if you have questions or if you would like to
%   use this code.

% Dealing with varargin and defaults
%
%code to allow for using a single cell array input
nargs=nargin;
if(nargs==1)
    varargin=varargin{1};
    nargs=length(varargin);
end

%single fluorophore intensity
sInt=1;

imSize = [512,512];
customIS=false;

%numerical aperature 100x lense
na = 1.49;

%wavelength of the three color channels
lambdaR = 675;%nm (RhodamineB emission maximum is 564)
lambdaG = 510;%nm (GFP emission maximum is 509)
lambdaB = 475;%nm

%points to use for fluorophore locations
pointsR = [];%vertcat(rand(10000,2)*511+1,[rand(40,1)*4+240,rand(40,1)*45+100]);
pointsG = [];
pointsB = [];

%Noise and standard deviation defaults
camNoise = 1;
camStdv = 0.1;
backNoise = 1;
backStdv = 0.1;

%simWindowSizeRatio
sdWindow = 5;

%pixelSize default
pixSize=64;

%illumination Mask
illuMask=ones(imSize);
customIM=false;

%overwrite defaults if name-value pair is provided
for p=1:2:nargs
    
    switch lower(varargin{p})
        case 'imsize'
            if numel(varargin{p+1}) == 2 && isreal(varargin{p+1}) && isinteger(varargin{p+1}) && all(varargin{p+1}>0)
                imSize=varargin{p+1};
                customIS=true;
            else
                error('The parameter "imSize" must be a two element array of positive real integers. Make sure the input is of type INT');
            end
            
        case 'lambdar'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>0
                lambdaR=varargin{p+1};
            else
                error('The parameter "lambdaR" must be a single positive real number');
            end
            
        case 'lambdag'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>0
                lambdaG=varargin{p+1};
            else
                error('The parameter "lambdaG" must be a single positive real number');
            end
            
        case 'lambdab'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>0
                lambdaB=varargin{p+1};
            else
                error('The parameter "lambdaB" must be a single positive real number');
            end
            
        case 'na'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>0
                na=varargin{p+1};
            else
                error('The parameter "na" must be a single positive real number');
            end
            
        case 'pointsr'
            if isreal(varargin{p+1})  && ( size(varargin{p+1},2)==2 || isempty(varargin{p+1}) )
                pointsR=varargin{p+1};
            else
               error('The parameter "pointsR" must be a nx2 array of real numbers, or an empty array'); 
            end
            
        case 'pointsg'
            if isreal(varargin{p+1})  && ( size(varargin{p+1},2)==2 || isempty(varargin{p+1}) )
                pointsG=varargin{p+1};
            else
               error('The parameter "pointsG" must be a nx2 array of real numbers, or an empty array'); 
            end
             
        case 'pointsb'
            if isreal(varargin{p+1})  && ( size(varargin{p+1},2)==2 || isempty(varargin{p+1}) )
                pointsB=varargin{p+1};
            else
               error('The parameter "pointsB" must be a nx2 array of real numbers, or an empty array'); 
            end
             
        case 'camnoise'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>=0
                camNoise=varargin{p+1};
            else
                error('The parameter "camNoise" must be a single real number, greater than or equal to zero');
            end
            
        case 'camstdv'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>=0
                camStdv=varargin{p+1};
            else
                error('The parameter "camStdv" must be a single real number, greater than or equal to zero');
            end
            
        case 'backnoise'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>=0
                backNoise=varargin{p+1};
            else
                error('The parameter "backNoise" must be a single real number, greater than or equal to zero');
            end
            
        case 'backstdv'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>=0
                backStdv=varargin{p+1};
            else
                error('The parameter "backStdv" must be a single real number, greater than or equal to zero');
            end
            
        case 'simwindowsizeratio'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>0
                sdWindow=varargin{p+1};
            else
                error('The parameter "simWindowSizeRatio" must be a single positive real number');
            end
            
        case 'pixsize'
            if numel(varargin{p+1}) == 1 && isreal(varargin{p+1}) && varargin{p+1}>0
                pixSize=varargin{p+1};
            else
                error('The parameter "pixSize" must be a single positive real number');
            end
        case 'illuminationmask'
            if numel(varargin{p+1}) >= 1 && all(isreal(varargin{p+1}))
                illuMask=varargin{p+1};
                customIM=true;
            else
                error('The parameter "illuminationMask" must be an array of real numbers.');
            end
        otherwise
            warning(['An unrecognize parameter was provided:',varargin{p}]);
    end
    
end

if(any(size(illuMask) ~= [imSize(2),imSize(1)]))
   illuMask= ones(imSize(2),imSize(1));
   if(customIS && customIM)
       warning('Your provided illumination mask"s size does not match your provided image size');
   elseif(customIM)
       error('Your provided illumination mask"s size does not match the default image size. Please provide a illumination mask of equal size to the default image size, or provide a custom image size.');
   else%customIM
       %They don't need to know about this one, it would be expected to
       %function correctly by creating the correct illuMaks
   end
end

%gaussian standard deviation (in pixels)
sigmaR=(lambdaR/(2*na))/pixSize;
sigmaG=(lambdaG/(2*na))/pixSize;
sigmaB=(lambdaB/(2*na))/pixSize;


outputImage = zeros(imSize(2),imSize(1),3);

%Camera Noise for all three channels
outputImage = outputImage + camNoise+camStdv.*randn(size(outputImage));

%Red
c=1;
%{
if(backNoise ~= 0 || backStdv ~= 0 )
    outputImage(:,:,c) = outputImage(:,:,c) + backNoise+backStdv.*randn(size(outputImage(:,:,c)));
end
%}
%tic;
for n=1:size(pointsR,1)
    rangeX = round(pointsR(n,1)-0.5*sdWindow*sigmaR):round(pointsR(n,1)+0.5*sdWindow*sigmaR);
    rangeX = rangeX(rangeX > 0 & rangeX <= imSize(1));
    rangeY = round(pointsR(n,2)-0.5*sdWindow*sigmaR):round(pointsR(n,2)+0.5*sdWindow*sigmaR);
    rangeY = rangeY(rangeY > 0 & rangeY <= imSize(2));
    if(numel(rangeX)>0 && numel(rangeY)>0)
        outputImage(rangeY,rangeX,c) = outputImage(rangeY,rangeX,c)+gaussIt2D(rangeX(1),rangeX(end),rangeY(1),rangeY(end),sigmaR,pointsR(n,1),pointsR(n,2));
    end
end
outputImage(:,:,c)=outputImage(:,:,c).*illuMask;
%{
display('Standard Gauss:');toc;
outputI2 = zeros(imSize(2),imSize(1),3);
tic;
outputI2(:,:,c) = outputI2(:,:,c) + gaussIt2DApprox(pointsR(:,1),pointsR(:,2),sigmaR,imSize,sdWindow*sigmaR,0.001);
display('Approx Gauss:');toc;
%}

%Green
c=2;
%{
if(backNoise ~= 0 || backStdv ~= 0 )
    outputImage(:,:,c) = outputImage(:,:,c) + backNoise+backStdv.*randn(size(outputImage(:,:,c)));
end
%}
for n=1:size(pointsG,1)
    rangeX = round(pointsG(n,1)-0.5*sdWindow*sigmaG):round(pointsG(n,1)+0.5*sdWindow*sigmaG);
    rangeX = rangeX(rangeX > 0 & rangeX <= imSize(1));
    rangeY = round(pointsG(n,2)-0.5*sdWindow*sigmaG):round(pointsG(n,2)+0.5*sdWindow*sigmaG);
    rangeY = rangeY(rangeY > 0 & rangeY <= imSize(2));
   if(numel(rangeX)>0 && numel(rangeY)>0)
        outputImage(rangeY,rangeX,c) = outputImage(rangeY,rangeX,c)+gaussIt2D(rangeX(1),rangeX(end),rangeY(1),rangeY(end),sigmaG,pointsG(n,1),pointsG(n,2));
    end
    %outputImage(rangeX,rangeY,c) = outputImage(rangeX,rangeY,c)+gaussIt2D(rangeX(1),rangeX(end),rangeY(1),rangeY(end),sigmaG,pointsG(n,1),pointsG(n,2));
end
outputImage(:,:,c)=outputImage(:,:,c).*illuMask;


%Blue
c=3;
%{
if(backNoise ~= 0 || backStdv ~= 0 )
    outputImage(:,:,c) = outputImage(:,:,c) + backNoise+backStdv.*randn(size(outputImage(:,:,c)));
end
%}
for n=1:size(pointsB,1)
    rangeX = round(pointsB(n,1)-0.5*sdWindow*sigmaB):round(pointsB(n,1)+0.5*sdWindow*sigmaB);
    rangeX = rangeX(rangeX > 0 & rangeX <= imSize(1));
    rangeY = round(pointsB(n,2)-0.5*sdWindow*sigmaB):round(pointsB(n,2)+0.5*sdWindow*sigmaB);
    rangeY = rangeY(rangeY > 0 & rangeY <= imSize(2));
    outputImage(rangeY,rangeX,c) = outputImage(rangeY,rangeX,c)+gaussIt2D(rangeX(1),rangeX(end),rangeY(1),rangeY(end),sigmaB,pointsB(n,1),pointsB(n,2));
end
outputImage(:,:,c)=outputImage(:,:,c).*illuMask;



end