%gaussIt2D
function [result]=gaussIt2D(xInit,xFinal,yInit,yFinal,sigma,xCenter,yCenter)
%Returns gaussian intenisty values for the region defined by 
%[xInit:xFinal,yInit:yFinal] using the gaussian properties sigma,centerX,centerY
[gridX,gridY]=meshgrid(xInit:xFinal,yInit:yFinal);

result=exp( -( (gridX-xCenter).^2 + (gridY-yCenter).^2 ) ./ (2*sigma.^2) );

end