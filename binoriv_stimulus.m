% Set Parameters corresponding to your environment and whatever you want
% Create a stimulus and save as an image file
% !Currently assuming left:right=red:blue
% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function [max_intensity_l, max_intensity_r, lineWidth, ann_rect, lineColour] = binoriv_stimulus()

% Uncomment out in case you want to run screen function, but your pc is not powerful to
% run with the error '----- ! PTB - ERROR: SYNCHRONIZATION FAILURE ! ----'
Screen('Preference', 'SkipSyncTests', 1);

mkdir stimulus

%% Parameters
savefile = 1; % 1: save data
orientation = [0 1 2 3]; % orientation: 0 for vertical, 1 for horizontal, 2 for binocular (you can maximise contrast for both gratings), 3 for overlapped image of 1&2
annulus = 1; % 0: no annulus; 1: annulus
annulus_save = 1;
photos = 0; % 0: no photos overlayed; 1: photos overlayed

% grating related
dist_scr = 42; % distance from screen (cm)
radius_gra = 5; % radius of grating circle (deg)
% grating's max intensity: defined based on maxmising blue's intensity
% (https://en.wikipedia.org/wiki/Relative_luminance)
max_intensity_l = 86.598;%255/2; % max contrast intensity [0 255]; 255 is the strongest
max_intensity_r = 255; 
bgint_l = 0; % min grating contrast intensity [0 255]
bgint_r = 0; 
cycles = 2; % spatial frequency (cycles per stimulus)
phase_l = -0.65; %5 % 20.5: FP on grating; 5: FP on background
phase_r = 20.5; %5 % 20.5: FP on grating; 5: FP on background


%% Define grating stimulus
prompt = 'How many inches is the screen?: ';
screen_inch = input(prompt);

GratDurationSecs = 2; % Abort demo after x seconds.
precue = 0.01*16.7*30; %16.7*30=501ms - 16.7ms is reflesh rate.


%%
AssertOpenGL;

% Get the list of screens and choose the one with the highest screen number.
screenNumber=max(Screen('Screens'));
%screenNumber=max(Screen('Screens')-1);

% Find the color values which correspond to white and black.
white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);


%% Define grating 
[w, rect] = Screen('OpenWindow', screenNumber, [0 0 0]);
HideCursor;
% make stimulus image
% if visual angle is less than 10°, tanV(deg) = S(size of
% stimulus)/D(distance from screen), i.e. S = D * tanV
D = dist_scr;   % viewing distance (cm)
V = radius_gra; % radius of grating circle (deg)
diameter = tand(V)*D; % cm
% from cm to px
%xysize = rect(4);
screen_diagonal = sqrt(rect(3)^2 + rect(4)^2);
xysize = diameter/(2.54/(screen_diagonal/screen_inch));
xylim = 2*pi*cycles;

% annulus
lineWidth = 15; % outline thickness
xysize_annulus = xysize/2 + lineWidth;
ann_rect = [rect(3)/2-xysize_annulus rect(4)/2-xysize_annulus rect(3)/2+xysize_annulus rect(4)/2+xysize_annulus];
lineColour = [255 255 255]; % outline colour

% colour of gratings
image_mono_v = zeros(round(xysize),round(xysize),3);
image_mono_h = zeros(round(xysize),round(xysize),3);
image_bino = zeros(round(xysize),round(xysize),3);
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
image_bino(:,:,1) = (bgint_l + (sin(x+phase_l)+1)*(max_intensity_l-bgint_l)/2) .* circle;
image_bino(:,:,3) = (bgint_r + (sin(y+phase_r)+1)*(max_intensity_r-bgint_r)/2) .* circle;

mono_v = Screen('MakeTexture', w, image_mono_v);
mono_h = Screen('MakeTexture', w, image_mono_h);
bino = Screen('MakeTexture', w, image_bino);

red_filter = zeros(rect(4),rect(3),3);
blue_filter = zeros(rect(4),rect(3),3);
red_filter(:,:,1) = 255;
red_filter(:,:,4) = 255*0.5;
red_filter = Screen('MakeTexture', w, red_filter);
blue_filter(:,:,3) = 255;
blue_filter(:,:,4) = 255*0.5;
blue_filter = Screen('MakeTexture', w, blue_filter);

%{
if savefile == 1 
    save('stimulus/grating_v.mat','image_mono_v')
    save('stimulus/grating_h.mat','image_mono_h')
end
%}

%% Load face and object images
%{
face = imread('stimulus/face.jpg');
objA = imread('stimulus/objectA.jpg');

face = imresize(face, [160 NaN]);
face(:,:,1) = face(:,:,1)*0.5;
face(:,:,2) = 0;
face(:,:,3) = 0;
objA = imresize(objA, [110 NaN]);
objA(:,:,1) = 0;
objA(:,:,2) = 0;
image_face = Screen('MakeTexture', w, face);
image_objA = Screen('MakeTexture', w, objA);
%}

%% Save images
[vbl, start] = Screen('Flip', w);

filename = [];
if annulus_save == 1
    Screen('FrameOval', w, lineColour , ann_rect, lineWidth);
    Screen('Flip', w); % update screen
    if savefile == 1
            imageArray = Screen('GetImage', w); % get colour information of the entire screen
            imwrite(imageArray, 'stimulus/annulus.png'); 
    end
end

for i = orientation
    if i == 0 % vertical
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('DrawTexture', w, mono_v)
        if photos == 1 
            Screen('DrawTexture', w, image_face)
        end
        %Screen('FillOval', w, [255/2 0 0], [1.5915*1000 0.8495*1000    1.6085*1000    0.8665*1000]);
        if annulus == 1
            Screen('FrameOval', w, lineColour , ann_rect, lineWidth);
        end
        %Screen('DrawTexture', w, red_filter)
        Screen('Flip', w); % update screen
        if savefile == 1
            imageArray = Screen('GetImage', w); % get colour information of the entire screen
            imwrite(imageArray, 'stimulus/left.png'); 
        end
    elseif i == 1 % horizontal
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('DrawTexture', w, mono_h)
        if photos == 1
            Screen('DrawTexture', w, image_objA)
        end
        %Screen('FillOval', w, [0 0 255], [1633.5 891.5 1650.5 908.5]);
        if annulus == 1
            Screen('FrameOval', w, lineColour , ann_rect, lineWidth);
        end
        %Screen('DrawTexture', w, blue_filter)
        Screen('Flip', w); % update screen
        if savefile == 1
            imageArray = Screen('GetImage', w); % get colour information of the entire screen
            imwrite(imageArray, 'stimulus/right.png'); 
        end
    elseif i == 2 % binocular
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('DrawTexture', w, bino)
        if photos == 1
            face(:,:,4) = 255;
            objA(:,:,4) = 255*0.5;
            image_face = Screen('MakeTexture', w, face);
            image_objA = Screen('MakeTexture', w, objA);
            Screen('DrawTexture', w, image_face)
            Screen('DrawTexture', w, image_objA)
        end
        if annulus == 1
            Screen('FrameOval', w, lineColour , ann_rect, lineWidth);
        end
        %Screen('DrawTexture', w, blue_filter)
        Screen('Flip', w); % update screen
        if savefile == 1
            imageArray = Screen('GetImage', w); % get colour information of the entire screen
            imwrite(imageArray, 'stimulus/binocular.png'); 
        end
    elseif i == 3 % overlapped image
        imdata_l  = imread('stimulus/left.png', 'png'); % vertical corresponds to left
        imdata_l(:,:,4) = 255; % alpha channel
        imdata_r  = imread('stimulus/right.png', 'png'); % horizontal corresponds to right
        imdata_r(:,:,4) = 255*0.5; % alpha channel
        imagetex_l = Screen('MakeTexture', w, imdata_l);
        imagetex_r = Screen('MakeTexture', w, imdata_r);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('DrawTexture', w, imagetex_l)
        Screen('DrawTexture', w, imagetex_r)
        Screen('Flip', w); % update screen
        if savefile == 1
            imageArray = Screen('GetImage', w); % get colour information of the entire screen
            imwrite(imageArray, 'stimulus/overlap.png'); 
        end
    end
end
 
Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
% The same commands wich close onscreen and offscreen windows also close textures.
sca;

fprintf('Well done! \n')
