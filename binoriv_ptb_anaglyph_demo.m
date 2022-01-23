% binoriv_ptb_anaglyph_demo
% (rotating) binocular rivalry stimulus
% Modified from http://keith.psych.udel.edu/ptb/part5.pdf
% Key presses start rotation, stop rotation and exit.
% Wear red/green glasses and fixate at the center.

cycles = 3; % spatial frequency (cycles per stimulus)
period = 1; % rotation period (sec)
max_red = 100;
max_blue = 230;
fix_red = 0;
fix_blue = 0;
fix_size = 5;

screenNumber = max(Screen('Screens'));
[w, rect] = Screen('OpenWindow', screenNumber, [0 0 0]);
HideCursor;
xc = rect(3)/2; % coordinates of screen center
yc = rect(4)/2;
% make stimulus image
xysize = rect(4)/4;
xylim = 2*pi*cycles;
[x,y] = meshgrid(-xylim:2*xylim/(xysize-1):xylim, ...
    -xylim:2*xylim/(xysize-1):xylim);
img = zeros(xysize,xysize,3);
circle = x.^2 + y.^2 <= xylim^2; % circular boundry
img(:,:,1) = 0+(sin(x)+1)/2*max_red .* circle; % red channel
% img(:,:,2) = (sin(y)+1)/2*255 .* circle; % green channel
img(:,:,3) = 0+(sin(y)+1)/2*max_blue .* circle; % blue channel (comment green and uncomment blue if googles are red/blue)

img(:,:,4) = 150;

t = Screen('MakeTexture', w, img);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('DrawTexture', w, t);
Screen('FillOval',w,[fix_red 0 0 255],[xc-fix_size-45 yc-fix_size-10 xc+fix_size-45 yc+fix_size-10]); %fixation pt. red
Screen('FillOval',w,[0 0 fix_blue 255],[xc-fix_size+45 yc-fix_size+15 xc+fix_size+45 yc+fix_size+15]); %fixation pt. blue



Screen('Flip', w);
KbWait;

% rotation part
%{
while KbCheck; end % make sure all keys are released
% animation loop
keyisdown = 0;
start_time = GetSecs;
while (~keyisdown)
    th = mod(360*(GetSecs-start_time)/period,360); % rotation angle
    Screen('DrawTexture', w, t, [], [], th);
    Screen('FillOval', w, [255 255 255], [xc-4 yc-4 xc+4 yc+4]);
    Screen('Flip', w);
    keyisdown = KbCheck;
end
while KbCheck; end
KbWait;
%}

ShowCursor;
Screen('Close',w);