function ax=subPlot(rows,columns,index,spacer,varargin)
%subPlot(rows,columns,index,spacer)
%
%Sets axis to subplots on a grid defined by 'rows' and 'columns'.
%'index' assigns the current axis to a subplot in position 'index' row-wise
%'spacer' is an optional border space between adjecent axis and between the
%axis and the edge of the figure.Default 'spacer' is 0.025 (2.5%) of the 
%figure size.
%
%
%If you want a cleaner display, remove axis ticks with the code below
%set(gca,'xtick',[]);set(gca,'ytick',[]);
%after plotting to the axis set up by subPlot
%
%Author: Taylor A. Reid (University of Minnesota)
%Date created: 10/14/2014
%Distributable under BSD liscence. 
%Copyright (c) 2014 Taylor A. Reid.
%All rights reserved.

extraTopSpacer = 0; 

if(nargin<4)
    spacer=0.025;
end
if(nargin>4)
    for ii = 1:2:numel(varargin)
        switch varargin{ii}
            case 'ExtraTopSpacer'
                extraTopSpacer = varargin{ii+1};
        end
    end
end

x=1;%total normalized x size
y=1;%total normalized y size
w=(x-(columns+1)*spacer)/columns;
h=(y-(rows+1)*(spacer+extraTopSpacer))/rows;

left=repmat(0+spacer*(1:columns)+w*(0:(columns-1)),[1,rows]);
bottom=repmat(0+spacer*(rows:-1:1)+h*((rows-1):-1:0),[columns,1]);

ax=subplot('Position',[left(index),bottom(index),w,h]);

end