%This is to take data from the EB1 dynamics model and then plot it to
%create a video showing the growth and then you can create a kymograph from
%it. 
 %addpath('Z:\cbs_lab_klei0091\Sam\Taylor code for MT dynamics with tpx2\MicroTubule\ImageSim\ImageSim');
%addpath('Z:\cbs_lab_klei0091\Sam\Taylor code for MT dynamics with tpx2\MicroTubule\Include');
 %addpath('R:\cbs_lab_klei0091\Sam\Taylor code for MT dynamics with tpx2\MicroTubule\ImageSim\ImageSim');
 %addpath('R:\cbs_lab_klei0091\Sam\Taylor code for MT dynamics with tpx2\MicroTubule\Include');
for i=1
   
    index=1;
    number=num2str(i);
    proteinLocations=readmatrix('proteinBindingAndRemovals1Condition1.xlsx','Sheet',i);
    protoLengths=readmatrix('protofilamentLengths1Condition1.xlsx','Sheet',i);
 videoOutputName=strcat('protofilamentLengths1Condition1',num2str(i),'-20FPSBackground1For64nMPixelsEvery5Frames.avi');%Give the name for the output video
 aviobj1 = VideoWriter(videoOutputName, 'Uncompressed AVI');
    aviobj1.FrameRate = 20;  %frame rate in frames/sec
   
    open(aviobj1);  %now opening the avi object
                        for j=1:5:size(protoLengths,1)%-1 is because protolengths has titles in row 1
                            proteins=proteinLocations(:,2*j-1:2*j);
                           outputImage=ImageSimControllerFunc(proteins,protoLengths(j+1,1:13));%plus 1 is cause protolenghs has titles in row 1
                            figure;
%         subPlot(1,1,1);
        imshow(mat2gray(outputImage)) ;
                            frame = getframe(gcf);   %captures CURRENT image
    writeVideo(aviobj1,frame);
                            close()
                        end
                        close(aviobj1);
                        index=index+1;
                        close all;
                        
end