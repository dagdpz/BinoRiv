% Set Parameters corresponding to your environment and whatever you want
% Create a stimulus and save as an image file
% Author: Ryo Segawa (whizznihil.kid@gmail.com)

% Uncomment out in case you want to run screen function, but your pc is not powerful to
% run with the error '----- ! PTB - ERROR: SYNCHRONIZATION FAILURE ! ----'
Screen('Preference', 'SkipSyncTests', 1);

mkdir stimulus

%% Parameters
savefile = 1; % 1: save data
orientation = [0 1]; % orientation: 0 for vertical, 1 for horizontal

dist_scr = 42; % distance from screen (cm)
radius_gra = 5; % radius of grating circle (deg)
cont = 255; % contrast intensity [1 255]; 255 is the strongest

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
cycles = 2; % spatial frequency (cycles per stimulus)
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
image_mono_v = zeros(round(xysize),round(xysize),3);
image_mono_h = zeros(round(xysize),round(xysize),3);
%img_r = zeros(xysize,xysize,3);
%img_g = zeros(xysize,xysize,3);
[x,y] = meshgrid(-xylim:2*xylim/(xysize-1):xylim, ...
    -xylim:2*xylim/(xysize-1):xylim);
circle = x.^2 + y.^2 <= xylim^2; % circular boundry
% create gratings
image_mono_v(:,:,1) = (sin(x)+1)/2*cont .* circle; % R channel; 255 is max contrast
%image_mono_v(:,:,2) = (sin(x)+1)/2*cont .* circle; % G channel
%image_mono_v(:,:,3) = (sin(x)+1)/2*cont .* circle; % B channel
image_mono_v(:,:,4) = 255/2; % alpha chanel (transparency); the drawing destination in perceptual
%image_mono_h(:,:,1) = (sin(y)+1)/2*cont .* circle;
%image_mono_h(:,:,2) = (sin(y)+1)/2*cont .* circle;
image_mono_h(:,:,3) = (sin(y)+1)/2*cont .* circle;
image_mono_h(:,:,4) = 255/2; % alpha chanel (transparency)
%image_mono_v = (sin(x)+1)/2*cont .* circle; % 256 is max contrast
%image_mono_h = (sin(y)+1)/2*cont .* circle;
%img_r(:,:,1) = (sin(x)+1)/2*255 .* circle; % red channel
%img_g(:,:,2) = (sin(y)+1)/2*255 .* circle; % green channel
% img(:,:,3) = (sin(y)+1)/2*255 .* circle; % blue channel (comment green and uncomment blue if googles are red/blue)
gratingtex_mono_v = Screen('MakeTexture', w, image_mono_v);
gratingtex_mono_h = Screen('MakeTexture', w, image_mono_h);
%gratingtex_r = Screen('MakeTexture', w, img_r);
%gratingtex_g= Screen('MakeTexture', w, img_g);

if savefile == 1 
    save('stimulus/grating_v.mat','image_mono_v')
    save('stimulus/grating_h.mat','image_mono_h')
end

%% Animation loop: Run until timeout
[vbl, start] = Screen('Flip', w);

filename = [];
for i = orientation
    if i == 0 % vertical
        Screen('DrawTexture', w, gratingtex_mono_v)
        %Screen('DrawTexture', w, gratingtex_r)
        Screen('Flip', w); % 画面を更新
        if savefile == 1
            imageArray = Screen('GetImage', w); % 画面全体の色情報を取得
            imwrite(imageArray, 'stimulus/grating_v.png'); 
        end
    elseif i == 1 % horizontal
        Screen('DrawTexture', w, gratingtex_mono_h)
%        Screen('DrawTexture', w, gratingtex_g)
        Screen('Flip', w); % 画面を更新
        if savefile == 1
            imageArray = Screen('GetImage', w); % 画面全体の色情報を取得
            imwrite(imageArray, 'stimulus/grating_h.png'); 
        end
    end
end

Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')

% The same commands wich close onscreen and offscreen windows also close textures.
sca;


fprintf('Well done! \n')
