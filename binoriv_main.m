% Set Parameters corresponding to your environment and whatever you want
% Author: Ryo Segawa (whizznihil.kid@gmail.com)

clear all
close all

% Uncomment out in case you want to run screen function, but your pc is not powerful to run with the error '----- ! PTB - ERROR: SYNCHRONIZATION FAILURE ! ----'
Screen('Preference', 'SkipSyncTests', 1);

AssertOpenGL;
screenNumber=max(Screen('Screens')); % use largest screen if using multiple displays
%screenNumber=max(Screen('Screens')-1);
[w, rect] = Screen('OpenWindow', screenNumber, [0 0 0]);
white=WhiteIndex(screenNumber); % Find the color values which correspond to white
red = [255 0 0];
blue = [0 0 255];
HideCursor;

%% Parameters
subj_type = 0; % 0 is human, 1 is monkey
report = 1; %0 is no report, 1 is report (i.e. record key pressing)

fix_size = 0.25; % diameter of a fixation spot (deg)
%colour_comb = 0; % 0 is (Red Blue), 1 is (Blue Red)
num_superblock = 5; % the number of super-blocks

screen_inch = 15.6; % size of the screen (inch)
dist_scr = 42; % distance from screen (cm)

%% Ask subject's name and num of attempts for output file
prompt = 'What is the subject''s first name?: ';
subj = input(prompt,'s');
if  ischar(subj)
else
    !echo Please type only characters.
    return
end

subj_dist = fullfile('recording/', subj);
mkdir(subj_dist, '/report/phys')
mkdir(subj_dist, '/report/bino')

%% Define fixation point
% if visual angle is less than 10°, tanV(deg) = S(size of stimulus)/D(distance from screen), i.e. S = D * tanV
D = dist_scr;   % viewing distance (cm)
V = fix_size; % radius of grating circle (deg)
diameter = tand(V)*D; % cm
% from cm to px
%xysize = rect(4);
screen_diagonal = sqrt(rect(3)^2 + rect(4)^2);
screen_inch = screen_inch;
xysize = diameter/(2.54/(screen_diagonal/screen_inch));
radius = xysize/2; % px

[centre(1), centre(2)] = RectCenter(rect);
theta = 1; % distance of a fixation spot from the centre (deg)
fix_d = 2 * dist_scr * tand(theta/2) * (180/pi); % distance of a fixation spot from the centre (rad)
fix_cent = [centre(1:1)-radius centre(2:2)-radius centre(1:1)+radius centre(2:2)+radius];
fix_left = [centre(1:1)-fix_d-radius centre(2:2)-radius centre(1:1)-fix_d+radius centre(2:2)+radius];
fix_right = [centre(1:1)+fix_d-radius centre(2:2)-radius centre(1:1)+fix_d+radius centre(2:2)+radius];
fix_below = [centre(1:1)-radius centre(2:2)-fix_d-radius centre(1:1)+radius centre(2:2)-fix_d+radius];
fix_up = [centre(1:1)-radius centre(2:2)+fix_d-radius centre(1:1)+radius centre(2:2)+fix_d+radius];
potential_loc = [fix_left; fix_right; fix_up; fix_below];

%% define other parameters
if subj_type == 0 % human
    trial_len = 2.0; % 2000 ms
    num_trial = 8;
    num_triad = 12;
    rest = 8; % break time (sec) after physical/binocular blocks % 8000 ms
    brk = 120; % break time (sec) after one super-block 
elseif subj_type == 1 % monkey
    trial_len = 0.8; % 800 ms
    num_trial = 8;
    num_triad = 12;
end
% shuffle the order of the phys stimuli
phys_stim = [];
for i=1:num_trial/2;phys_stim=horzcat(phys_stim,0);end
for i=1:num_trial/2;phys_stim=horzcat(phys_stim,1);end

%% load stimuli
%{
imdata_l = load('stimulus/grating_v','image_mono_v');
imdata_l = imdata_l.image_mono_v;
imdata_r = load('stimulus/grating_h','image_mono_h');
imdata_r = imdata_r.image_mono_h;
%}
imdata_l  = imread('stimulus/left.png', 'png'); % vertical corresponds to left
imdata_l(:,:,4) = 255*0.5; % alpha channel
imdata_r  = imread('stimulus/right.png', 'png'); % horizontal corresponds to right
imdata_r(:,:,4) = 255*0.5; % alpha channel
imagetex_l = Screen('MakeTexture', w, imdata_l);
imagetex_r = Screen('MakeTexture', w, imdata_r);

%% Main
for i = 1:num_superblock
            
    for j = 1:num_triad
        phys_bino = randi([0 1],1,1); % if 0 phys->bino, if 1 bino->phys
        phys_stim = phys_stim(randperm(num_trial)); % shuffle the order of the phys stimuli
        if phys_bino == 0  
            phys_switch(w,red,blue,i,j,trial_len,num_trial,imagetex_l,imagetex_r,potential_loc,report,phys_stim,subj_dist);
            percep_switch(w,red,blue,i,j,trial_len,num_trial,imagetex_l,imagetex_r,potential_loc,report,subj_dist);
        elseif phys_bino == 1
            percep_switch(w,red,blue,i,j,trial_len,num_trial,imagetex_l,imagetex_r,potential_loc,report,subj_dist);
            phys_switch(w,red,blue,i,j,trial_len,num_trial,imagetex_l,imagetex_r,potential_loc,report,phys_stim,subj_dist);
        end

        fprintf('Rest now... \n')
        [vbl, start] = Screen('Flip', w);
        while (vbl < (rest + start))
            vbl = Screen('Flip', w); % return current time
        end
        Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
    end
    
    % Break
    fprintf('Break now... \n')
    [vbl, start] = Screen('Flip', w);
    while (vbl < (brk + start))
        vbl = Screen('Flip', w); % return current time
    end
    
end