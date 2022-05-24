% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function binoriv_phys_switch(w,superblock,triad)

global VAR
global SETTINGS

if VAR.eye_track == 1
    durTrialMax = VAR.trial_len; % max trial duration [s]
    fsMatlab = 1000; % matlab "sampling" frequency [Hz]
    bufferSize = durTrialMax * fsMatlab;
    memoryBuffer_fixloc = []; % [fix_loc]
    memoryBuffer_eyepos = []; % buffer data: [time trial x_eye y_eye]
    memoryBuffer_eyesize = []; % buffer data: [time trial width hight]
    counter = 0;
end
    
if VAR.report == 1
    % Define the key for reporting the content of perception
    KbName('UnifyKeyNames');
    VertKey = KbName('UpArrow');
    HorzKey = KbName('LeftArrow');
    
    vert_press = [];
    horz_press = [];
    pre_VertKey = 0;
    pre_HorzKey = 0;
end

% fix_loc_label = randi([1 4],1,VAR.num_trial);
fix_loc_label = randi([1 2],1,VAR.num_trial);

answer = [];

%% animation
[vbl, start] = Screen('Flip', w);
% for online plotting
if VAR.eye_track == 1; figure('Name','Online GUI','Position',SETTINGS.GUI_coordinates,'Renderer', 'Painters'); axis([-5 5 -5 5]); axis square; grid on; hold on; end 
for t = 1:VAR.num_trial
    if VAR.grating_task == 1
        if VAR.phys_stim(t) == 1 % Right
%             fix_loc = VAR.potential_loc_grat(fix_loc_label(1,t),:);
            fix_loc = VAR.potential_loc_grat(4,:);
        else %Left
%             fix_loc = VAR.potential_loc(fix_loc_label(1,t),:);
            fix_loc = VAR.potential_loc(1,:);
        end
    else
        fix_loc = VAR.potential_loc(fix_loc_label(1,t),:);
    end
    
    
    if VAR.phys_stim(t) == 0 % L
        answer = vertcat(answer, 'L');
    elseif VAR.phys_stim(t) == 1 % R
        answer = vertcat(answer, 'R');
    end

    while (vbl < (t*VAR.trial_len + start))
        if VAR.phys_stim(t) == 0 % L
            if VAR.grating_task == 1
                Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 0 0 0]);
                Screen('DrawTexture', w, VAR.imagetex_l);
                Screen('FillOval', w, VAR.red, fix_loc);
                Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 1]);
                Screen('FrameOval', w, VAR.lineColour , VAR.ann_rect, VAR.lineWidth); % annulus
            else
                Screen('FillOval', w, VAR.red, fix_loc);
            end
        elseif VAR.phys_stim(t) == 1 % R
            if VAR.grating_task == 1
                Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [0 0 1 0]);
                Screen('DrawTexture', w, VAR.imagetex_r);
                Screen('FillOval', w, VAR.blue, fix_loc);
                Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 1]);
                Screen('FrameOval', w, VAR.lineColour , VAR.ann_rect, VAR.lineWidth); % annulus
            else
                Screen('FillOval', w, VAR.blue, fix_loc);
            end
        end
        vbl = Screen('Flip', w); % return current time
        
        [keyIsDown, press, KeyCode] = KbCheck;
        if KeyCode(VAR.EscapeKey)==1 
            break
        end
        if VAR.report == 1
            % Record the time if subject pressed button
            if KeyCode(VertKey)==1 && pre_VertKey==0
                pre_VertKey = 1;
                vert_press = vertcat(vert_press, (vbl-start)*1000); % ms
            elseif KeyCode(VertKey)==0 && pre_VertKey==1
                pre_VertKey = 0;
                vert_press = vertcat(vert_press, (vbl-start)*1000); % ms
            elseif KeyCode(HorzKey)==1 && pre_HorzKey==0
                pre_HorzKey = 1;
                horz_press = vertcat(horz_press, (vbl-start)*1000); % ms
            elseif KeyCode(HorzKey)==0 && pre_HorzKey==1
                pre_HorzKey = 0;
                horz_press = vertcat(horz_press, (vbl-start)*1000); % ms
            end
        end
        if VAR.eye_track == 1
            [x_eye, y_eye] = aux_GetCalibratedEyeHandPos();
            [width,hight] = vpx_GetPupilSize();
            counter = counter + 1;
            %memoryBuffer_fixloc(counter,1:4) = fix_loc;
            memoryBuffer_fixloc(counter,1) = fix_loc_label(1,t);
            memoryBuffer_eyepos(counter,1) = (vbl-start)*1000;
            memoryBuffer_eyepos(counter,2) = t;
            memoryBuffer_eyepos(counter,3) = x_eye;
            memoryBuffer_eyepos(counter,4) = y_eye;
            memoryBuffer_eyesize(counter,1) = (vbl-start)*1000;
            memoryBuffer_eyesize(counter,2) = t;
            memoryBuffer_eyesize(counter,3) = width;
            memoryBuffer_eyesize(counter,4) = hight;
            %memoryBuffer_eyepos(counter,:) = [(vbl-start)*1000,t,x_eye,y_eye];
            %memoryBuffer_eyesize(counter,:) = [(vbl-start)*1000,t,width,hight];
            
            % Online plot
            try 
                x_eye_pre = memoryBuffer_eyepos(counter-1,3);
                y_eye_pre = memoryBuffer_eyepos(counter-1,4);
                line([x_eye_pre x_eye], [y_eye_pre y_eye], 'Color','k','Marker','none','LineWidth',0.2)
            catch
            end
            drawnow
        end
        if VAR.eye_track == 1
            filename = [VAR.fig_dir '/reco_online_phys_' num2str(superblock) '_' num2str(triad) '_' num2str(t) '.png'];
            saveas(gcf,filename)
            filename = [VAR.fig_dir '/reco_online_phys_' num2str(superblock) '_' num2str(triad) '_' num2str(t) '.fig'];
            saveas(gcf,filename)
        end
    end
end

%% save
if VAR.report == 1
    vert_press = array2table(vert_press);
    horz_press = array2table(horz_press);
    filename = [VAR.subj_dist '/report/phys/vertpress_repo_' num2str(superblock) '_' num2str(triad) '.csv'];
    writetable(vert_press, filename, 'WriteVariableNames', false);
    filename = [VAR.subj_dist '/report/phys/horzpress_repo_' num2str(superblock) '_' num2str(triad) '.csv'];
    writetable(horz_press, filename, 'WriteVariableNames', false);
    answer = array2table(answer);
    filename = [VAR.subj_dist '/report/phys/answer_' num2str(superblock) '_' num2str(triad) '.csv'];
    writetable(answer, filename, 'WriteVariableNames', false);
end



if VAR.eye_track == 1
    if VAR.report == 1
        filename = [VAR.subj_dist '/report/phys/fixloc_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_fixloc);
        writetable(buffer, filename, 'WriteVariableNames', false);
        filename = [VAR.subj_dist '/report/phys/eyepos_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_eyepos);
        writetable(buffer, filename, 'WriteVariableNames', false);
        filename = [VAR.subj_dist '/report/phys/eyesize_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_eyesize);
        writetable(buffer, filename, 'WriteVariableNames', false);
    else
        filename = [VAR.norepo_dir '/phys/fixloc_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_fixloc);
        writetable(buffer, filename, 'WriteVariableNames', false);
        filename = [VAR.norepo_dir '/phys/eyepos_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_eyepos);
        writetable(buffer, filename, 'WriteVariableNames', false);
        filename = [VAR.norepo_dir '/phys/eyesize_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_eyesize);
        writetable(buffer, filename, 'WriteVariableNames', false);
    end
end