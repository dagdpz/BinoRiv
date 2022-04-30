% Set Parameters corresponding to your environment and whatever you want
% Create a stimulus and save as an image file
% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function [max_intensity_l, max_intensity_r, lineWidth, ann_rect, lineColour] = binoriv_stimulus()

global VAR

% Uncomment out in case you want to run screen function, but your pc is not powerful to
% run with the error '----- ! PTB - ERROR: SYNCHRONIZATION FAILURE ! ----'
Screen('Preference', 'SkipSyncTests', 1);

mkdir stimulus

%% Parameters
screen_inch = 27; % screen size (inch)
screen_width = 61-1.4;
screen_hight = 33.7;
savefile = 1; % 1: save data
orientation = [0 1 2 3 4]; % orientation: 0 for vertical, 1 for horizontal, 2 for binocular (you can maximise contrast for both gratings), 3 for overlapped image of 1&2
annulus = 1; % 0: no annulus; 1: annulus
annulus_save = 1;
photos = 0; % 0: no photos overlayed; 1: photos overlayed
subj = 'Luba_0412';

% grating-related
dist_scr = 47; % distance from screen (cm)
diameter_gra = 10%5; % diameter of grating circle (deg)
sf = 0.5%2; % spatial frequency (cycle/degree)
%num_cycle = 10; % spatial frequency (cycles/stimulus)% grating's max intensity: defined based on maxmising blue's intensity
% (https://en.wikipedia.org/wiki/Relative_luminance)
path_luminance_l = ['stimulus/' subj '/luminance_l_' subj '.mat'];
try path_luminance_r = ['stimulus/' subj '/luminance_r_' subj '.mat']; catch; end
if exist(path_luminance_l, 'file') == 2 && exist(path_luminance_r, 'file') == 2 
    max_intensity_l = load(path_luminance_l);
    max_intensity_r = load(path_luminance_r);
    max_intensity_l = 150%max_intensity_l.max_intensity_l;
    max_intensity_r = 255%max_intensity_r.max_intensity_r;
else
    max_blue_intensity = 150;
    max_intensity_l = 0.3396*max_blue_intensity; % max contrast intensity [0 255]; 255 is the strongest
    max_intensity_r = max_blue_intensity; 
end
%gtline_width_l = 0%0.8*max_intensity_l;
%gtline_width_r = 0%0.8*max_intensity_r;
bgint_l = 0;%0.7*max_intensity_l; % min grating contrast intensity [0 255]
bgint_r = 0;%0.7*max_intensity_r; 
phase = (1/2)*pi;
phase_l = phase;%20.5; % 20.5: FP on grating; 4.5: FP on background for cycle=2
phase_r = phase;%20.5%4.6 %20.5: FP on grating; 4.5: FP on background for cycle=2


%%
AssertOpenGL;

% Get the list of screens and choose the one with the highest screen number.
screenNumber=max(Screen('Screens'));
%screenNumber=max(Screen('Screens')-1);


%% Define grating 
[w, rect] = Screen('OpenWindow', screenNumber, [0 0 0]);
HideCursor;
% make stimulus image
% if visual angle is less than 10°, tanV(deg) = S(size of stimulus)/D(distance from screen), i.e. S = D * tanV
D = dist_scr;   % viewing distance (cm)
V = diameter_gra; % diameter of grating circle (deg)
diameter = 2*D * tand(V/2);%tand(V)*D % cm
% from cm to px:
VAR.width_pxpercm = rect(3)/screen_width;
hight_cxpercm = rect(4)/screen_hight;
xysize = diameter * VAR.width_pxpercm; % px
%adj = 2*radius_gra/num_cycle; % in case the unit of spatial frequency is cycles/stimulus
cycles = sf * diameter_gra; % in case the unit of spatial frequency is cycles/degree
%cycles = radius_gra/adj;%3.8; % in case the unit of spatial frequency is cycles/stimulus
xylim = pi*cycles;

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
%[x,y] = meshgrid(-xylim:2*xylim/(rect(4)-1):xylim, ...
 %   -xylim:2*xylim/(rect(4)-1):xylim);
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


%{
image_bino_ind(:,:,1) = (sin(x_ind)+1)/2*max_intensity_l .* circle; 
image_bino_ind(:,:,3) = (sin(y_ind)+1)/2*max_intensity_r .* circle;
image_mask(:,:,:) = 255;
%}

mono_v = Screen('MakeTexture', w, image_mono_v);
mono_h = Screen('MakeTexture', w, image_mono_h);
bino = Screen('MakeTexture', w, image_bino);
%bino = Screen('MakeTexture', w, image_bino_ind);
%mask = Screen('MakeTexture', w, image_mask);

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


%% pick up coordinates at max value of gratings 
VAR.wave_length = (xysize/(cycles));
VAR.wave_length_half = VAR.wave_length/2;
VAR.x_cent = rect(3)/2;
VAR.y_cent = rect(4)/2;

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
        %Screen('FillOval', w, [51.9588 0 0], 1.0e+03 *[1.5702    0.8967    1.5769    0.9033]);
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
        %Screen('FillOval', w, [0 0 153], 1.0e+03 *[1.5887    0.8547    1.5953    0.8613]);
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
        
        %{
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 0 0 0]);
        Screen('FillOval', w, [255 0 0], [1229.48551330144	716.486706422681	1236.51210045608	723.513293577319]);
        Screen('FillOval', w, [255 0 0], [1323.48789954392	716.486706422681	1330.51448669856	723.513293577319]);
        Screen('FillOval', w, [255 0 0], [1276.48670642268	763.487899543919	1283.51329357732	770.514486698558]);
        Screen('FillOval', w, [255 0 0], [1276.48670642268	669.485513301442	1283.51329357732	676.512100456081]);
        %
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [0 0 1 0]);
        Screen('FillOval', w, [0 0 255], [924.665900601040	551.165900601040	928.728240338313	555.228240338313]);
        Screen('FillOval', w, [0 0 255], [1006.66590060104	660.165900601040	1010.72824033831	664.228240338313]);
        Screen('FillOval', w, [0 0 255], [924.665900601040	660.165900601040	928.728240338313	664.228240338313]);
        %Screen('FillOval', w, [0 0 255], [1006.66590060104	551.165900601040	1010.72824033831	555.228240338313]);
        %}
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
        
    elseif i == 4 % debug
        imdata_l(:,:,4) = 255; % alpha channel
        imdata_r(:,:,4) = 255; % alpha channel
        imagetex_l = Screen('MakeTexture', w, imdata_l);
        imagetex_r = Screen('MakeTexture', w, imdata_r);
        %
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 0 0 0]);
        Screen('DrawTexture', w, imagetex_l)
        %Screen('FillOval', w, [0 0 0], [1229.48551330144	716.486706422681	1236.51210045608	723.513293577319]);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [0 0 1 0]);
        Screen('DrawTexture', w, imagetex_r)
        %Screen('FillOval', w, [0 0 0], [1276.48670642268	669.485513301442	1283.51329357732	676.512100456081]);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 0]);
%        FP = max_loc(1,3)
        Screen('FillOval', w, [255 255 255], [rect(3)/2-5 rect(4)/2-5 rect(3)/2+5 rect(4)/2+5]);
        wave_length = (xysize/(cycles*2));
        Screen('FillOval', w, [255 255 255], [rect(3)/2-wave_length-3 rect(4)/2-wave_length-3 rect(3)/2-wave_length+3 rect(4)/2-wave_length+3]);
        %Screen('FillOval', w, [255 255 255], [VAR.image_left_edge+max_loc_index-3 VAR.image_top_edge+max_loc_index-3 VAR.image_left_edge+max_loc_index+3 VAR.image_top_edge+max_loc_index+3]);
        %Screen('FillOval', w, [255 0 0], fix_upleft_l);
        %Screen('FillOval', w, [0 0 255], fix_upleft_r);
        Screen('Flip', w); % update screen
        if savefile == 1
            imageArray = Screen('GetImage', w); % get colour information of the entire screen
            imwrite(imageArray, 'stimulus/debug.png'); 
        end
    end
end
 

Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
% The same commands wich close onscreen and offscreen windows also close textures.
sca;

fprintf('Stimuli created! \n')
