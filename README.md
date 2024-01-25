# GenerateSimulatedMTImages
Generate simulated microtubules with EB1 tip tracking from the EB1 simulation. 

## About the project
This code is used to generate simulated microtubule images with EB1 tip tracking that are generated from the EB1 Dynamics Simulation script. Of note, this code will take as input the matched proteinBindingAndRemovals... and protofilamentLengths... files that are generated as output from the EB1 Dynamics Simulation Script. Then, this script will read in the information about the protofilaments and the EB1 positions to generate a simulated EB1 image for each timestep recorded from the simulation. This will be saved as an output AVI file that can be used for downstream analysis (such as making kymographs or measuring tip specificity). 


### Built With
MATLAB by Mathworks

## Getting Started

To use the code, download the script into the folder with your excels shees for proteinBindingAndRemovals and protofilamentLengths that were output from your dynamic EB1 simulation. Of note, each sheet of these files will be a different simulation run; this script will generate videos for all of the sheets if you adjust the parameters for the initial i for loop to include all sheets. Next, if desired, you can adjust the associated pixel size (in the ImageSimControllerFunc). By default, the pixel size is 64 nm which is consisent with our 100X objective on our TIRF microscope. The output will be an avi movie with the simulated images. 

## Prerequisites

Ensure Matlab is operational on your device. 

## Installation

This script requires some add-ons including: 
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox


## Usage

The way the code works is that it labels every single EB1 with a green fluorophore. It then labels a small percentage of the tubulin dimers with a fluorophore (similar to how in these types of experimental assays, not all tubulin dimers are labelled with a fluorphore). It will then plot these points and use a PSF to generate an image similar to what would be seen from a TIRF microscope. Of note, noise is added to the image to more closely mimic what is seen experimentally. 


Of note, this code was used for some publications from the Gardner lab including: 

https://www.biorxiv.org/content/10.1101/2022.06.07.495114v1.abstract 




## Contact

Samuel Gonzalez-(https://www.linkedin.com/in/samuel-gonzalez-081504163/) - samueljgonzalez@hotmail.com
