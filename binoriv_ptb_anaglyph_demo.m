% (rotating) binocular rivalry stimulus
% Modified from http://keith.psych.udel.edu/ptb/part5.pdf
% Key presses start rotation, stop rotation and exit.
% Wear red/green glasses and fixate at the center.

cycles = 20; % spatial frequency (cycles per stimulus)
period = 1; % rotation period (sec)
screenNumber = max(Screen('Screens'));
[w, rect] = Screen('OpenWindow', screenNumber, [0 0 0]);
HideCursor;
xc = rect(3)/2; % coordinates of screen center
yc = rect(4)/2;
% make stimulus image
xysize = rect(4);
xylim = 2*pi*cycles;
[x,y] = meshgrid(-xylim:2*xylim/(xysize-1):xylim, ...
    -xylim:2*xylim/(xysize-1):xylim);
img = zeros(xysize,xysize,3);
circle = x.^2 + y.^2 <= xylim^2; % circular boundry
img(:,:,1) = (sin(x)+1)/2*255 .* circle; % red channel
img(:,:,2) = (sin(y)+1)/2*255 .* circle; % green channel
% img(:,:,3) = (sin(y)+1)/2*255 .* circle; % blue channel (comment green and uncomment blue if googles are red/blue)

t = Screen('MakeTexture', w, img);
Screen('DrawTexture', w, t);
Screen('FillOval',w,[255 255 255],[xc-4 yc-4 xc+4 yc+4]); %fixation pt.
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