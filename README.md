# GenerateSimulatedMTImages
Generate simulated microtubules with EB1 tip tracking from the EB1 simulation. 

## About the project
This code is used to generate simulated microtubule images with EB1 tip tracking that are generated from the EB1 Dynamics Simulation script. Of note, this code will take as input the matched proteinBindingAndRemovals... and protofilamentLengths... files that are generated as output from the EB1 Dynamics Simulation Script. Then, this script will read in the information about the protofilaments and the EB1 positions to generate a simulated EB1 image for each timestep recorded from the simulation. This will be saved as an output AVI file that can be used for downstream analysis (such as making kymographs or measuring tip specificity). Of note, there are two main files to used as the top script file labelled as Gonzalez_simulated_MT_growth and Gonzalez_simulated_MT_growthForEvL. Essentially, their inputs are different. If you ran the original batching code for simulated EB1 tip tracking with the modality to track binding initially at edge versus lattice sites, then you should use the script with EvL. Otherwise, you need to use the script ending in _growth. The reason for this is that the input is slightly different and the difference in the scripts accounts for that. Of note, for EvL, the green channel will show EB1 at edges while the blue channel will show EB1 at lattice sites; it was done this way so they could easily be seperated in the output image but the blue channel is a bit hard to see. 


### Built With
MATLAB by Mathworks

## Getting Started
To get started, download this Github repository as a zipped file by downloading it from the "code" dropdown menu, unzip the file locally, and save it in an accessible location.  

To use the code, move the repository files into the folder with your excels sheets for proteinBindingAndRemovals and protofilamentLengths that were output from your dynamic EB1 simulation. Of note, each sheet of these files will be a different simulation run; this script will generate videos for all of the sheets if you adjust the parameters for the initial i for loop to include all sheets. Next, if desired, you can adjust the associated pixel size (in the ImageSimControllerFunc or ImageSimControllerFuncEvL depending on which modality you are using). By default, the pixel size is 64 nm which is consisent with our 100X objective on our TIRF microscope. The output will be an avi movie with the simulated movies that can be used for downstream analysis. 

## Prerequisites

Ensure Matlab is operational on your device. This code was generated to work with MATLAB 2022b. It should work with newer versions of Matlab. 

## Installation

This script requires some add-ons including: 
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox


## Usage

The way the code works is that it labels every single EB1 with a green fluorophore. It then labels a small percentage of the tubulin dimers with a fluorophore (similar to how in these types of experimental assays, not all tubulin dimers are labelled with a fluorphore). It will then plot these points and use a PSF to generate an image similar to what would be seen from a TIRF microscope. Of note, noise is added to the image to more closely mimic what is seen experimentally. To get started, download all matlab files into the same folder as your excel files that are output from the SimulatedEB1TipTracking repository (specifically, the files that start with proteinbBindingAndRemovals... and protofilamentLengths.... Next, open up the Gonzalez_simulated_MT_Growth function (or the one ending in EvL if you want to have the EBs that bound initially to edge sites to be differentially colored than those that initially bound to lattice sites). Next, in this file (Gonzalez_simulated_MT_Growth...), make sure to change the names of the excel files that are being loaded as input, choose the output file name, and adjust the frames per rate for the video output if you so desire (standard is 20 frames per minute). Finally, also adjust the main for loop (i) so that you can decide how many of the simulation runs will be generated into movies (each simulation run will be kept in a different sheet within the two input files). The code will then generate simulated images for every  5th recorded timesteps from the input files (this can gain be adjusted if so desired and tends to be about every 5-7 seconds of  run time in the simulation)


Of note, this code was used for some publications from the Gardner lab including: 

https://www.biorxiv.org/content/10.1101/2022.06.07.495114v1.abstract 




## Contact

Samuel Gonzalez-(https://www.linkedin.com/in/samuel-gonzalez-081504163/) - samueljgonzalez@hotmail.com
