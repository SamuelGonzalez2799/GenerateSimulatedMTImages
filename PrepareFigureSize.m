% Cleanup and prepare for maximazation
function answer = PrepareFigureSize

set(0,'Units','pixels')  % this and the next 2 lines set up the program to "almost" maximize images to the screensize
a = get(0,'ScreenSize');
set(0,'DefaultFigurePosition',[a(3)*.1,a(4)*.1,a(3)*.8,a(4)*.8]);

answer = 1;

end