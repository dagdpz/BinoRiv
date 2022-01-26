% Set Parameters corresponding to your environment and whatever you want
% Author: Ryo Segawa (whizznihil.kid@gmail.com)

clear all
close all

% delete later
[max_intensity_l, max_intensity_r, lineWidth, ann_rect, lineColour] = binoriv_stimulus();

% Uncomment out in case you want to run screen function, but your pc is not powerful to run with the error '----- ! PTB - ERROR: SYNCHRONIZATION FAILURE ! ----'
Screen('Preference', 'SkipSyncTests', 1);

AssertOpenGL;
screenNumber=max(Screen('Screens')); % use largest screen if using multiple displays
%screenNumber=max(Screen('Screens')-1);

[w, rect] = Screen('OpenWindow', screenNumber, [0 0 0]);
white=WhiteIndex(screenNumber); % Find the color values which correspond to white

HideCursor;

KbName('UnifyKeyNames');
EscapeKey = KbName('ESCAPE');

%% Parameters
subj = 'Ryo'; % subject's name
subj_type = 0; % 0 is human, 1 is monkey
report = 0; %0 is no report (record only eye-tracking), 1 is report (i.e. record also key pressing)

fix_size = 0.15%0.25; % diameter of a fixation spot (deg)
theta = 1.0; % distance of a fixation spot from the centre (deg)
%colour_comb = 0; % 0 is (Red Blue), 1 is (Blue Red)
% Contrast of the fixed point to the maximum luminance of the grating; luminance contrast is defined as Weber contrast (https://en.m.wikipedia.org/wiki/Contrast_(vision))
% contrast = (intensity of FP - intensity of background)/ intensity of background 
contrast = -0.5; % [-1.0 1.0]; e.g. -0.5 means luminance of 50% to the maximum luminance of the grating

screen_inch = 24; % size of the screen (inch)
dist_scr = 42; % distance from screen (cm)

%% duration-related parameters
num_superblock = 1%5; % the number of super-blocks
if subj_type == 0 % human
    trial_len = 2.0; % 2000 ms
    num_trial = 20%8;
    num_triad = 6%12;
    rest = 0%8; % break time (sec) after physical/binocular blocks % 8000 ms
    brk = 0%120; % break time (sec) after one super-block 
    subj_dist = fullfile('recording/human', subj);
elseif subj_type == 1 % monkey
    trial_len = 0.8; % 800 ms
    num_trial = 8;
    num_triad = 12;
    subj_dist = fullfile('recording/monkey', subj);
end

% shuffle the order of the phys stimuli
phys_stim = [];
for i=1:num_trial/2;phys_stim=horzcat(phys_stim,0);end
for i=1:num_trial/2;phys_stim=horzcat(phys_stim,1);end

%% Create folders
if report == 1
    mkdir(subj_dist, '/report/phys')
    mkdir(subj_dist, '/report/bino')
end

%% Define fixation point
% intensity of the fixation points
red_intensity = luminance(contrast, "red", max_intensity_l, max_intensity_r);
red = [red_intensity 0 0]; 
%red = [51.9588 0 0];
blue_intensity = luminance(contrast, "blue", max_intensity_l, max_intensity_r);
blue = [0 0 blue_intensity];
%blue = [0 0 153];

% if visual angle is less than 10°, tanV(deg) = S(size of stimulus)/D(distance from screen), i.e. S = D * tanV
D = dist_scr;   % viewing distance (cm)
V = fix_size; % radius of grating circle (deg)
diameter = tand(V)*D; % cm
% from cm to px
%xysize = rect(4);
screen_diagonal = sqrt(rect(3)^2 + rect(4)^2);
screen_inch = screen_inch;
xysize = diameter/(2.54/(screen_diagonal/screen_inch)); % 2.54 is cm/inch
radius = xysize/2; % px

[centre(1), centre(2)] = RectCenter(rect);
fix_d = 2 * dist_scr * tand(theta/2) * (180/pi); % distance of a fixation spot from the centre (rad)
fix_cent = [centre(1:1)-radius centre(2:2)-radius centre(1:1)+radius centre(2:2)+radius];
fix_left = [centre(1:1)-fix_d-radius centre(2:2)-radius centre(1:1)-fix_d+radius centre(2:2)+radius];
fix_right = [centre(1:1)+fix_d-radius centre(2:2)-radius centre(1:1)+fix_d+radius centre(2:2)+radius];
fix_below = [centre(1:1)-radius centre(2:2)-fix_d-radius centre(1:1)+radius centre(2:2)-fix_d+radius];
fix_up = [centre(1:1)-radius centre(2:2)+fix_d-radius centre(1:1)+radius centre(2:2)+fix_d+radius];
potential_loc = [fix_left; fix_right; fix_up; fix_below];

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
%imdata_l(:,:,4) = 255*0.5; % alpha channel % image should be lower contrast than FP!
imdata_r  = imread('stimulus/right.png', 'png'); % horizontal corresponds to right
%imdata_r(:,:,4) = 255*0.5; % alpha channel
imagetex_l = Screen('MakeTexture', w, imdata_l);
imagetex_r = Screen('MakeTexture', w, imdata_r);
imdata_bino  = imread('stimulus/binocular.png', 'png'); % binocular image 
%imdata_bino(:,:,4) = 255*0.5; % alpha channel
imagetex_bino = Screen('MakeTexture', w, imdata_bino);

imdata_annulus  = imread('stimulus/annulus.png', 'png'); % annulus
imagetex_annulus = Screen('MakeTexture', w, imdata_annulus);


%% Main
for i = 1:num_superblock
    
    for j = 1:num_triad
        phys_bino = 1%randi([0 1],1,1); % if 0 phys->bino, if 1 bino->phys
        phys_stim = phys_stim(randperm(num_trial)); % shuffle the order of the phys stimuli
        if phys_bino == 0
            binoriv_phys_switch(w,red,blue,i,j,trial_len,num_trial, ...
                imagetex_l,imagetex_r, ...
                potential_loc,report,phys_stim, ...
                subj_dist,EscapeKey, ...
                lineWidth, ann_rect, lineColour);
            binoriv_percep_switch(w,red,blue,i,j,trial_len,num_trial, ...
                imagetex_l,imagetex_r,imagetex_bino, ...
                potential_loc,report, ...
                subj_dist,EscapeKey, ...
                lineWidth, ann_rect, lineColour);
        elseif phys_bino == 1
            binoriv_percep_switch(w,red,blue,i,j,trial_len,num_trial, ...
                imagetex_l,imagetex_r,imagetex_bino, ...
                potential_loc,report, ...
                subj_dist,EscapeKey, ...
                lineWidth, ann_rect, lineColour);
            %binoriv_phys_switch(w,red,blue,i,j,trial_len,num_trial, ...
             %   imagetex_l,imagetex_r, ...
              %  potential_loc,report,phys_stim, ...
               % subj_dist,EscapeKey, ...
                %lineWidth, ann_rect, lineColour);
        end
        
        fprintf('Rest now... \n')
        [vbl, start] = Screen('Flip', w);
        while (vbl < (rest + start))
            vbl = Screen('Flip', w); % update screen
        end
        Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
    end
    
    % Break
    fprintf('Break now... \n')
    [vbl, start] = Screen('Flip', w);
    while (vbl < (brk + start))
        [keyIsDown, press, KeyCode] = KbCheck;
        if KeyCode(EscapeKey)==1 
            break
        end
        vbl = Screen('Flip', w); % return current time
    end
    
end

ShowCursor;
Screen('Close',w);


function [I_fx] = luminance(contrast, colour, max_intensity_l, max_intensity_r)
    if colour == "red"
        I_gr = max_intensity_l;
    elseif colour == "blue"
        I_gr = max_intensity_r;
    end
    I_fx = I_gr*(2+contrast)/(2-contrast);
end