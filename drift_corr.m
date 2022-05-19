% Fix gaze offsets automatically
% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function drift_corr(w, rect, radius)
global VAR

screen_width = 61-1.4;
dist_scr = 47; % cm

[vbl, start] = Screen('Flip', w);
x_eye_list =[];
y_eye_list =[];
while (vbl-start) < 2.0
    Screen('FillOval', w, [255 255 255], [rect(3)/2-radius rect(4)/2-radius rect(3)/2+radius rect(4)/2+radius]);
    vbl = Screen('Flip', w); % return current time
    
    [x_eye, y_eye] = aux_GetCalibratedEyeHandPos(); % deg
    x_eye_array = vertcat(x_eye_list, x_eye);
    y_eye_array = vertcat(y_eye_list, y_eye);
    
    % from px to deg|fixation spot
%     x_eye_cent = x_eye - rect(3)/2; % px from centre
%     y_eye_cent = y_eye - rect(4)/2;
%     width_cmperpx = screen_width/rect(3); % cm/px of the screen
%     x_eye_cent = x_eye_cent*width_cmperpx; % cm from centre
%     y_eye_cent = y_eye_cent*width_cmperpx;
%     x_eye_cent = 2*atand(x_eye_cent/(2*dist_scr)); % deg
%     y_eye_cent = 2*atand(y_eye_cent/(2*dist_scr));
%     distance = sqrt(x_eye_cent^2+y_eye_cent^2);

    [keyIsDown, press, KeyCode] = KbCheck;
    if KeyCode(VAR.EscapeKey)==1 
        break
    end
end

x_eye_mean = mean(x_eye_array);
y_eye_mean = mean(y_eye_array);

if x_eye_mean > 0.25 || x_eye_mean < 0.25
    VAR.eye_offset_x = VAR.eye_offset_x - x_eye_mean; % deg
elseif y_eye_mean > 0.25 || y_eye_mean < 0.25
    VAR.eye_offset_y = VAR.eye_offset_y - y_eye_mean; % deg
end