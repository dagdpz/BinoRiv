# BinoRiv
This repository provides the binocular rivalry task with anaglyph used in the project ''Probing the Link between Perceptual Visibility and Neural Correlates of Consciousness''.  
Contact: Ryo Segawa (rsegawa@dpz.eu/whizznihil.kid@gmail.com)

## Prerequisites
* MATLAB (R2018b or later version)
* Psychophysics Toolbox Version 3 (PTB-3): http://psychtoolbox.org/
* Viewpoint eye-tracking system
  
## Procedure
1. check_luminance_l.m, check_luminance_r.m  (for human)  
Upper and down arrows are corresponding to increase and decrease luminance. After adjusting the luminance, press esc key - then the luminance data is saved.  
2. Calibration by some ways, and save the offsets and gains  
3. binoriv_main.m  
Execute the task.  
Subfunction: binoriv_stimulus.m, phys_switch.m, percep_switch.m  
4. binoriv_repo_plot_linking_ver2.m  
Plot results of the task - stimulus switches, report outcome, gaze position.  
Subfunction: square_colouring.m (depends on later ver of Matlab)  
4. raw_eye_pos.m
Output same results as binoriv_repo_plot_linking_ver2.m + raw gaze position in a square diagram  
Subfunction: square_colouring.m (depends on later ver of Matlab)  
