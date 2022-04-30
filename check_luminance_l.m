% Create a stimulus, operate luminance
% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function [max_intensity_l, max_intensity_r, lineWidth, ann_rect, lineColour] = check_luminance_l()

% Uncomment out in case you want to run screen function, but your pc is not powerful to
% run with the error '----- ! PTB - ERROR: SYNCHRONIZATION FAILURE ! ----'
Screen('Preference', 'SkipSyncTests', 1);

% Define the key for reporting the content of perception
KbName('UnifyKeyNames');
UpKey = KbName('UpArrow');
DownKey = KbName('DownArrow');
EscapeKey = KbName('ESCAPE');


%% Parameters
screen_inch = 27; % screen size (inch)
annulus = 1; % 0: no annulus; 1: annulus
subj = 'Luba_0412';

% grating-related
dist_scr = 47; % distance from screen (cm)
radius_gra = 5; % radius of grating circle (deg)
sf = 2; % spatial frequency (cycle/degree)
%num_cycle = 10; % spatial frequency (cycles/stimulus)% grating's max intensity: defined based on maxmising blue's intensity
% (https://en.wikipedia.org/wiki/Relative_luminance)
max_blue_intensity = 150;
max_intensity_l = 0.3396*max_blue_intensity; % max contrast intensity [0 255]; 255 is the strongest
max_intensity_r = max_blue_intensity; 
%gtline_width_l = 0%0.8*max_intensity_l;
%gtline_width_r = 0%0.8*max_intensity_r;
bgint_l = 0;%.3*max_intensity_l; % min grating contrast intensity [0 255]
bgint_r = 0;%0.3*max_intensity_r; 
phase = (1/2)*pi;
phase_l = phase;%20.5; % 20.5: FP on grating; 4.5: FP on background for cycle=2
phase_r = phase;%20.5%4.6 %20.5: FP on grating; 4.5: FP on background for cycle=2


%%
AssertOpenGL;

% Get the list of screens and choose the one with the highest screen number.
screenNumber=max(Screen('Screens'));
%screenNumber=max(Screen('Screens')-1);

mkdir stimulus
mkdir(['stimulus/' subj])

%% Define grating 
[w, rect] = Screen('OpenWindow', screenNumber, [0 0 0]);
HideCursor;
% make stimulus image
% if visual angle is less than 10°, tanV(deg) = S(size of
% stimulus)/D(distance from screen), i.e. S = D * tanV
D = dist_scr;   % viewing distance (cm)
V = radius_gra; % radius of grating circle (deg)
diameter = tand(V)*D; % cm
% from cm to px:
dots = sqrt(rect(3)^2 + rect(4)^2);
inch = 27
dpi = dots/inch % px
xysize_true = diameter * dpi/2.54 % 2.54 is cm/inch
screen_diagonal = sqrt(rect(3)^2 + rect(4)^2);
xysize = diameter/(2.54/(screen_diagonal/screen_inch));
%adj = 2*radius_gra/num_cycle; % in case the unit of spatial frequency is cycles/stimulus
cycles = sf * radius_gra; % in case the unit of spatial frequency is cycles/degree
%cycles = radius_gra/adj;%3.8; % in case the unit of spatial frequency is cycles/stimulus
xylim = 2*pi*cycles;

% annulus
lineWidth = 15; % outline thickness
xysize_annulus = xysize/2 + lineWidth;
ann_rect = [rect(3)/2-xysize_annulus rect(4)/2-xysize_annulus rect(3)/2+xysize_annulus rect(4)/2+xysize_annulus];
lineColour = [255 255 255]; % outline colour
                
% colour of gratings
image_mono_v = zeros(round(xysize),round(xysize),3);
image_mono_h = zeros(round(xysize),round(xysize),3);
image_bino = zeros(round(xysize),round(xysize),3); % spatial freq depends on the size of stimulus
image_bino_ind = zeros(rect(4),rect(4),3); % spatial freq is NOT dependent on the size of stimulus
image_mask = zeros(rect(4),rect(4),3);
[x,y] = meshgrid(-xylim:2*xylim/(round(xysize)-1):xylim, ...
    -xylim:2*xylim/(round(xysize)-1):xylim);
circle = x.^2 + y.^2 <= xylim^2; % circular boundry
% create gratings
%image_mono_v(:,:,1) = 255; % R channel; 255 is max contrast
image_mono_v(:,:,1) =  (bgint_l + (sin(x+phase_l)+1)*(max_intensity_l-bgint_l)/2) .* circle; %(sin(x)+1)/2*max_intensity_l .* circle; % R channel; 255 is max contrast
%image_mono_v(:,:,2) = (sin(x)+1)/2*cont_l .* circle; % G channel
%image_mono_v(:,:,3) = (sin(x)+1)/2*cont_l .* circle; % B channel
%image_mono_v(:,:,4) = 255*0.5; % alpha chanel (transparency); the drawing destination in perceptual
%image_mono_h(:,:,1) = (sin(y)+1)/2*cont_r .* circle;
%image_mono_h(:,:,2) = (sin(y)+1)/2*cont_r .* circle;
%image_mono_h(:,:,3) = 255;
image_mono_h(:,:,3) = (bgint_r + (sin(y+phase_r)+1)*(max_intensity_r-bgint_r)/2) .* circle; %(sin(y)+1)/2*max_intensity_r .* circle;
%image_mono_h(:,:,4) = 255*0.5; % alpha chanel (transparency)
%image_bino(:,:,1) = (bgint_l + (sin(x+phase_l)+1)*(max_intensity_l-bgint_l-gtline_width_l)/2+gtline_width_l) .* circle;
%image_bino(:,:,3) = (bgint_r + (sin(y+phase_r)+1)*(max_intensity_r-bgint_r-gtline_width_r)/2+gtline_width_r) .* circle;
image_bino(:,:,1) = (bgint_l + (sin(x+phase_l)+1)*(max_intensity_l-bgint_l)/2) .* circle;
image_bino(:,:,3) = (bgint_r + (sin(y+phase_r)+1)*(max_intensity_r-bgint_r)/2) .* circle;


%% main
Screen('Flip', w);
[keyIsDown, press, KeyCode] = KbCheck;
while KeyCode(EscapeKey)==0
    image_mono_v(:,:,1) =  (bgint_l + (sin(x+phase_l)+1)*(max_intensity_l-bgint_l)/2) .* circle;
    mono_v = Screen('MakeTexture', w, image_mono_v);

    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('DrawTexture', w, mono_v)

    if annulus == 1
        Screen('FrameOval', w, lineColour , ann_rect, lineWidth);
    end
    Screen('Flip', w); % update screen

    [keyIsDown, press, KeyCode] = KbCheck;
    if KeyCode(UpKey)==1 && max_intensity_l < 255
        max_intensity_l = max_intensity_l + 1
    elseif KeyCode(DownKey)==1 && max_intensity_l >0
        max_intensity_l = max_intensity_l - 1
    end
end

filename = ['stimulus/' subj '/luminance_l_' subj '.mat'];
save(filename, 'max_intensity_l')
 
Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
% The same commands wich close onscreen and offscreen windows also close textures.
sca;
