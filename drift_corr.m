% Fix gaze offsets automatically
% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function drift_corr(w, superblock, triad)
global VAR
global SETTINGS
close all

[vbl, start] = Screen('Flip', w);
x_eye_array =[];
y_eye_array =[];

% for online plotting
% figure('Name','Online GUI','Position',SETTINGS.GUI_coordinates,'Renderer', 'Painters'); axis([-5 5 -5 5]); axis square; grid on; hold on; 
% plot(0,0,'*', 'MarkerEdgeColor','g')
counter = 0;
counter_corr = 0;

while (vbl-start) < VAR.dc_time
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 1]);
    Screen('FillOval', w, [255 255 255], [VAR.rect(3)/2-VAR.radius VAR.rect(4)/2-VAR.radius VAR.rect(3)/2+VAR.radius VAR.rect(4)/2+VAR.radius]);
    vbl = Screen('Flip', w); % return current time
    
    [x_eye, y_eye] = aux_GetCalibratedEyeHandPos(); % deg
    x_eye_array = vertcat(x_eye_array, x_eye);
    y_eye_array = vertcat(y_eye_array, y_eye);
    
    % from px to deg|fixation spot
%     x_eye_cent = x_eye - VAR.rect(3)/2; % px from centre
%     y_eye_cent = y_eye - VAR.rect(4)/2;
%     width_cmperpx = SETTINGS.screen_w_cm/VAR.rect(3); % cm/px of the screen
%     x_eye_cent = x_eye_cent*width_cmperpx; % cm from centre
%     y_eye_cent = y_eye_cent*width_cmperpx;
%     x_eye_cent = 2*atand(x_eye_cent/(2*SETTINGS.vd)); % deg
%     y_eye_cent = 2*atand(y_eye_cent/(2*SETTINGS.vd));
%     distance = sqrt(x_eye_cent^2+y_eye_cent^2);

    % Online plot
    counter = counter + 1;
%     try 
%         x_eye_pre = x_eye_array(counter-1,1);
%         y_eye_pre = y_eye_array(counter-1,1);
%         plot(x_eye,y_eye,'o', 'MarkerEdgeColor',[0.5 0.5 0.5])
%         line([x_eye_pre x_eye], [y_eye_pre y_eye], 'Color','k','Marker','none','LineWidth',0.2)
%     catch
%     end
%     drawnow
    
    
    % drift correction
    if (vbl-start) > VAR.dc_time/2 && (vbl-start) < (VAR.dc_time/2+VAR.dc_time/8) && counter_corr == 0 
        VAR.eye_offset_x = VAR.eye_offset_x - x_eye; % deg
        VAR.eye_offset_y = VAR.eye_offset_y - y_eye; % deg
        counter_corr = counter_corr + 1;
    elseif (vbl-start) > (VAR.dc_time/2+VAR.dc_time/8) && (vbl-start) < (VAR.dc_time/2+2*VAR.dc_time/8) && counter_corr == 1
        VAR.eye_offset_x = VAR.eye_offset_x - x_eye; % deg
        VAR.eye_offset_y = VAR.eye_offset_y - y_eye; % deg
        counter_corr = counter_corr + 1;
    elseif (vbl-start) > (VAR.dc_time/2+2*VAR.dc_time/8) && (vbl-start) < (VAR.dc_time/2+3*VAR.dc_time/8) && counter_corr == 2
        VAR.eye_offset_x = VAR.eye_offset_x - x_eye; % deg
        VAR.eye_offset_y = VAR.eye_offset_y - y_eye; % deg
        counter_corr = counter_corr + 1;
    elseif (vbl-start) > (VAR.dc_time/2+3*VAR.dc_time/8) && (vbl-start) < (VAR.dc_time/2+4*VAR.dc_time/8) && counter_corr == 3
        VAR.eye_offset_x = VAR.eye_offset_x - x_eye; % deg
        VAR.eye_offset_y = VAR.eye_offset_y - y_eye; % deg
        counter_corr = counter_corr + 1;
    end

    [keyIsDown, press, KeyCode] = KbCheck;
    if KeyCode(VAR.EscapeKey)==1 
        break
    end
end

% filename = [VAR.fig_dir '/reco_online_DC_' num2str(superblock) '_' num2str(triad) '.png'];
% saveas(gcf,filename)
% filename = [VAR.fig_dir '/reco_online_DC_' num2str(superblock) '_' num2str(triad) '.fig'];
% saveas(gcf,filename)

% x_eye_mean = mean(x_eye_array);
% y_eye_mean = mean(y_eye_array);
% 
% if x_eye_mean > 0.25 || x_eye_mean < 0.25
%     VAR.eye_offset_x = VAR.eye_offset_x - x_eye_mean; % deg
% elseif y_eye_mean > 0.25 || y_eye_mean < 0.25
%     VAR.eye_offset_y = VAR.eye_offset_y - y_eye_mean; % deg
% end