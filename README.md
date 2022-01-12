# BinoRiv
This repository provides the binocular rivalry task for anaglyph used in the project ''Probing the Link between Perceptual Visibility and Neural Correlates of Consciousness''.  
Contact: Ryo Segawa (rsegawa@dpz.eu/whizznihil.kid@gmail.com)

## Prerequisites
* MATLAB (R2018b or later version)
* Psychophysics Toolbox Version 3 (PTB-3): http://psychtoolbox.org/

## Concept
The code framework consists of a programme to create stimuli, a programme to read the created stimuli and present the task, and a programme to plot the results. The main execution flow of the task looks like this:
![code_str](https://user-images.githubusercontent.com/41120302/148266421-a6cf3ded-f9a8-45e3-a45d-64c135d4ccff.png)

## Procedure
1. binoriv_stimulus.m  
Create stimulus images for the task.
2. main.m  
Show the stimulus image created in the previous step and execute the task.  
Subfunction: phys_switch.m, percep_switch.m
3. binoriv_repo_plot.m/  
Plots the results of the task.  
Subfunction: square_colouring.m  
  
The project is ongoing and significant changes are to be expected...  

## Stimuli
![stimulus left](/stimulus/left.png)
![stimulus right](/stimulus/right.png)
![stimulus overlap](/stimulus/overlap.png)