% Set Parameters corresponding to your environment and whatever you want
% Author: Ryo Segawa (whizznihil.kid@gmail.com)

clear all
close all

global VAR
global PATH
global SETTINGS

% delete later
VAR.subj = 'Ryo_0524_DC'; % subject's name
[max_intensity_l, max_intensity_r, VAR.lineWidth, VAR.ann_rect, VAR.lineColour] = binoriv_stimulus();

% Uncomment out in case you want to run screen function, but your pc is not powerful to run with the error '----- ! PTB - ERROR: SYNCHRONIZATION FAILURE ! ----'
%Screen('Preference', 'SkipSyncTests', 1);

AssertOpenGL;
%screenNumber=max(Screen('Screens')); % use largest screen if using multiple displays
screenNumber=max(Screen('Screens')-1);

[w, rect] = Screen('OpenWindow', screenNumber, [0 0 0]);
white=WhiteIndex(screenNumber); % Find the color values which correspond to white

HideCursor;

KbName('UnifyKeyNames');
VAR.EscapeKey = KbName('ESCAPE');

%% Set the monitor-related parameters
SETTINGS.GUI_coordinates = [-1550 20 700 700*rect(4)/rect(3)];
SETTINGS.screen_h_cm = 33.7; % screen hight (cm)
SETTINGS.screen_w_cm = 61-1.4; % screen width (cm)
SETTINGS.vd = 47; % distance from screen (cm)
SETTINGS.screen_uh_cm = SETTINGS.screen_h_cm/2;
SETTINGS.screen_lh_deg = atan((SETTINGS.screen_h_cm - SETTINGS.screen_uh_cm)/SETTINGS.vd)/(pi/180);
SETTINGS.screen_uh_deg = atan(SETTINGS.screen_uh_cm/SETTINGS.vd)/(pi/180);
SETTINGS.screen_w_deg = 2*atan((SETTINGS.screen_w_cm/2)/SETTINGS.vd)/(pi/180);
SETTINGS.screen_h_deg       = SETTINGS.screen_lh_deg + SETTINGS.screen_uh_deg;
SETTINGS.matlab_version = datevec(version('-date'));

% Online GUI-related
SETTINGS.GUI_coordinates = [-2000 200 2000 2000*rect(4)/rect(3)];

%% Parameters
VAR.subj = 'Ryo_0524_DC'; % subject's name
VAR.subj_type = 0; % 0 is human, 1 is monkey
VAR.report = 1; % 0 is no report, 1 is report (i.e. record key pressing)
VAR.eye_track = 1; % 0: eye tracker on, 1: eye tracker off
VAR.grating_task = 1; % 0 if real task, 1 if background is only grating 

fix_size = 0.2; % diameter of a fixation spot (deg)
FP_loc = 1; % 0: fixed distance (theta) from the centre, 1: FP onto non-grating-overlapped place, 2: FP onto grating-overlapped place
distfromcent = 2; % integer
theta = 1.0; % distance of a fixation spot from the centre (deg)
colour_comb = 1; % 0 is (left:Red right:Blue), 1 is (left:Blue right:Red)
% Contrast of the fixed point to the maximum luminance of the grating; luminance contrast is defined as Weber contrast (https://en.m.wikipedia.org/wiki/Contrast_(vision))
% contrast = (intensity of FP - intensity of background)/ intensity of background 
contrast = -0.5; % [-1.0 1.0]; e.g. -0.5 means luminance of 50% to the maximum luminance of the grating

screen_inch = 27; % size of the screen (inch)
dist_scr = SETTINGS.vd; % distance from screen (cm)

PATH.calibpath = 'D:\Data\human'; % path where you've got calibration data

%% duration-related parameters
num_superblock = 1%5; % the number of super-blocks
if VAR.subj_type == 0 % human
    VAR.trial_len = 2.0; % 2000 ms
    VAR.num_trial = 8;
    num_triad = 12;
    rest = 0%8; % break time (sec) after physical/binocular blocks % 8000 ms
    brk = 0%120; % break time (sec) after one super-block 
    VAR.subj_dist = fullfile('recording/human', VAR.subj);
elseif VAR.subj_type == 1 % monkey
    VAR.trial_len = 0.8; % 800 ms
    VAR.num_trial = 8;
    num_triad = 12;
    rest = 8; % break time (sec) after physical/binocular blocks % 8000 ms
    brk = 120; % break time (sec) after one super-block 
    VAR.subj_dist = fullfile('recording/monkey', VAR.subj);
end

% shuffle the order of the phys stimuli
VAR.phys_stim = [];
for i=1:VAR.num_trial/2;VAR.phys_stim=horzcat(VAR.phys_stim,0);end
for i=1:VAR.num_trial/2;VAR.phys_stim=horzcat(VAR.phys_stim,1);end

%% Create folders
if VAR.report == 0
    mkdir(VAR.subj_dist, '/noreport/phys')
    mkdir(VAR.subj_dist, '/noreport/bino')
    VAR.norepo_dir = [VAR.subj_dist '/noreport'];
    VAR.fig_dir = [VAR.norepo_dir '/figures'];
    mkdir(VAR.fig_dir)
elseif VAR.report == 1    
    mkdir(VAR.subj_dist, '/report/phys')
    mkdir(VAR.subj_dist, '/report/bino')
    VAR.repo_dir = [VAR.subj_dist '/report'];
    VAR.fig_dir = [VAR.repo_dir '/figures'];
    mkdir(VAR.fig_dir)
end


%% Define fixation point
% intensity of the fixation points
red_intensity = luminance(contrast, 'red', max_intensity_l, max_intensity_r);
VAR.red = [red_intensity 0 0]; 
blue_intensity = luminance(contrast, 'blue', max_intensity_l, max_intensity_r);
VAR.blue = [0 0 blue_intensity];

% if visual angle is less than 10°, tanV(deg) = S(size of stimulus)/D(distance from screen), i.e. S = D * tanV
D = dist_scr;   % viewing distance (cm)
V = fix_size; % radius of grating circle (deg)
diameter = tand(V)*D; % cm
% from cm to px
%xysize = rect(4);
screen_diagonal = sqrt(rect(3)^2 + rect(4)^2);
screen_inch = screen_inch;
xysize = diameter * VAR.width_pxpercm;%diameter/(2.54/(screen_diagonal/screen_inch)); % 2.54 is cm/inch
radius = xysize/2; % px

% fixation point locations
wave_length_dist = VAR.wave_length*distfromcent;
if FP_loc == 1 %FP onto non-grating-overlapped place
    fix_upleft_l = [VAR.x_cent-wave_length_dist-radius VAR.y_cent-wave_length_dist-VAR.wave_length_half-radius VAR.x_cent-wave_length_dist+radius VAR.y_cent-wave_length_dist-VAR.wave_length_half+radius];
    fix_upright_l = [VAR.x_cent+wave_length_dist-radius VAR.y_cent-wave_length_dist-VAR.wave_length_half-radius VAR.x_cent+wave_length_dist+radius VAR.y_cent-wave_length_dist-VAR.wave_length_half+radius];
    fix_belowleft_l = [VAR.x_cent-wave_length_dist-radius VAR.y_cent+wave_length_dist+VAR.wave_length_half-radius VAR.x_cent-wave_length_dist+radius VAR.y_cent+wave_length_dist+VAR.wave_length_half+radius];
    fix_belowright_l = [VAR.x_cent+wave_length_dist-radius VAR.y_cent+wave_length_dist+VAR.wave_length_half-radius VAR.x_cent+wave_length_dist+radius VAR.y_cent+wave_length_dist+VAR.wave_length_half+radius];
    VAR.potential_loc = [fix_upleft_l; fix_upright_l; fix_belowleft_l; fix_belowright_l];
    
    fix_upleft_r = [VAR.x_cent-wave_length_dist-VAR.wave_length_half-radius VAR.y_cent-wave_length_dist-radius VAR.x_cent-wave_length_dist-VAR.wave_length_half+radius VAR.y_cent-wave_length_dist+radius];
    fix_upright_r = [VAR.x_cent+wave_length_dist+VAR.wave_length_half-radius VAR.y_cent-wave_length_dist-radius VAR.x_cent+wave_length_dist+VAR.wave_length_half+radius VAR.y_cent-wave_length_dist+radius];
    fix_belowleft_r = [VAR.x_cent-wave_length_dist-VAR.wave_length_half-radius VAR.y_cent+wave_length_dist-radius VAR.x_cent-wave_length_dist-VAR.wave_length_half+radius VAR.y_cent+wave_length_dist+radius];
    fix_belowright_r = [VAR.x_cent+wave_length_dist+VAR.wave_length_half-radius VAR.y_cent+wave_length_dist-radius VAR.x_cent+wave_length_dist+VAR.wave_length_half+radius VAR.y_cent+wave_length_dist+radius];
    VAR.potential_loc_grat = [fix_upleft_r; fix_upright_r; fix_belowleft_r; fix_belowright_r];

elseif FP_loc == 2 %FP onto grating-overlapped place
    fix_left = [VAR.x_cent-wave_length_dist-radius VAR.y_cent-wave_length_dist-radius VAR.x_cent-wave_length_dist+radius VAR.y_cent-wave_length_dist+radius];
    fix_up = [VAR.x_cent+wave_length_dist-radius VAR.y_cent-wave_length_dist-radius VAR.x_cent+wave_length_dist+radius VAR.y_cent-wave_length_dist+radius];
    fix_below = [VAR.x_cent-wave_length_dist-radius VAR.y_cent+wave_length_dist-radius VAR.x_cent-wave_length_dist+radius VAR.y_cent+wave_length_dist+radius];
    fix_right = [VAR.x_cent+wave_length_dist-radius VAR.y_cent+wave_length_dist-radius VAR.x_cent+wave_length_dist+radius VAR.y_cent+wave_length_dist+radius];
    VAR.potential_loc = [fix_left; fix_right; fix_up; fix_below];
    VAR.potential_loc_grat = [fix_left; fix_right; fix_up; fix_below];
else
    [centre(1), centre(2)] = RectCenter(rect);
    fix_d = 2 * dist_scr * tand(theta/2) * (180/pi); % distance of a fixation spot from the centre (left)
    fix_cent = [centre(1:1)-radius centre(2:2)-radius centre(1:1)+radius centre(2:2)+radius];

    fix_left = [centre(1:1)-fix_d-radius centre(2:2)-radius centre(1:1)-fix_d+radius centre(2:2)+radius];
    fix_right = [centre(1:1)+fix_d-radius centre(2:2)-radius centre(1:1)+fix_d+radius centre(2:2)+radius];
    fix_below = [centre(1:1)-radius centre(2:2)-fix_d-radius centre(1:1)+radius centre(2:2)-fix_d+radius];
    fix_up = [centre(1:1)-radius centre(2:2)+fix_d-radius centre(1:1)+radius centre(2:2)+fix_d+radius];
    VAR.potential_loc = [fix_left; fix_right; fix_up; fix_below];

    fix_d = fix_d/sqrt(2); % distance of a fixation spot from the centre (right)
    fix_upleft = [centre(1:1)-radius-fix_d centre(2:2)-radius+fix_d centre(1:1)+radius-fix_d centre(2:2)+radius+fix_d];
    fix_upright = [centre(1:1)-radius+fix_d centre(2:2)-radius+fix_d centre(1:1)+radius+fix_d centre(2:2)+radius+fix_d];
    fix_belowleft = [centre(1:1)-radius-fix_d centre(2:2)-radius-fix_d centre(1:1)+radius-fix_d centre(2:2)+radius-fix_d];
    fix_belowright = [centre(1:1)-radius+fix_d centre(2:2)-radius-fix_d centre(1:1)+radius+fix_d centre(2:2)+radius-fix_d];
    VAR.potential_loc_grat = [fix_upleft; fix_belowright; fix_upright; fix_belowleft];
end

% shuffle the order of the phys stimuli
VAR.phys_stim = [];
for i=1:VAR.num_trial/2;VAR.phys_stim=horzcat(VAR.phys_stim,0);end
for i=1:VAR.num_trial/2;VAR.phys_stim=horzcat(VAR.phys_stim,1);end

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
VAR.imagetex_l = Screen('MakeTexture', w, imdata_l);
VAR.imagetex_r = Screen('MakeTexture', w, imdata_r);
imdata_bino  = imread('stimulus/binocular.png', 'png'); % binocular image 
%imdata_bino(:,:,4) = 255*0.5; % alpha channel
VAR.imagetex_bino = Screen('MakeTexture', w, imdata_bino);

imdata_annulus  = imread('stimulus/annulus.png', 'png'); % annulus
imagetex_annulus = Screen('MakeTexture', w, imdata_annulus);

%% load calibration data
if VAR.eye_track == 1
    if exist([PATH.calibpath filesep 'last_eyecal.mat'],'file'),
        disp(['Found previous eye calibation in ' PATH.calibpath filesep 'last_eyecal.mat']);
        load([PATH.calibpath filesep 'last_eyecal.mat']);
        VAR.eye_gain_x = eye_gain_x;
        VAR.eye_gain_y = eye_gain_y;
        VAR.eye_offset_x = eye_offset_x;
        VAR.eye_offset_y = eye_offset_y;
    else
        fprintf('Have not been calibrated! \n')
        pause
    end
end

%% Main
for i = 1:num_superblock
    switch_log = [];
    for j = 1:num_triad
        phys_bino = 1%randi([0 1],1,1); % if 0 phys->bino, if 1 bino->phys
        VAR.phys_stim = VAR.phys_stim(randperm(VAR.num_trial)); % shuffle the order of the phys stimuli
       
        % Drift correction
        if VAR.eye_track == 1; drift_corr(w, rect, radius, i, j); end
        % Switches
        if phys_bino == 0
            binoriv_phys_switch(w,i,j);
            binoriv_percep_switch(w,i,j);
        elseif phys_bino == 1
            binoriv_percep_switch(w,i,j);
            binoriv_phys_switch(w,i,j);
        end
        switch_log = vertcat(switch_log, phys_bino);

        fprintf('Rest now... \n')
        [vbl, start] = Screen('Flip', w);
        %if VAR.eye_track == 1 || VAR.report == 1 ; try rtplot(i,j); catch; end; end
        
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
        if KeyCode(VAR.EscapeKey)==1 
            break
        end
        vbl = Screen('Flip', w); % return current time
    end
    
    % save switch log
    if VAR.report == 1
        filename = [VAR.subj_dist '/report/switch_' num2str(i) '.csv'];
    elseif VAR.report == 0
        filename = [VAR.subj_dist '/noreport/switch_' num2str(i) '.csv'];
    end
    switch_log = array2table(switch_log);
    writetable(switch_log, filename, 'WriteVariableNames', false)
    
end

ShowCursor;
Screen('Close',w);


%% save the variables
dt = datetime('today');
DateString = datestr(dt,'yyyymmdd');
if VAR.report == 0
    filename = [VAR.subj_dist '/noreport/variables_repo_' DateString '.mat'];
elseif VAR.report == 1
    filename = [VAR.subj_dist '/report/variables_repo_' DateString '.mat'];
end
save(filename)

