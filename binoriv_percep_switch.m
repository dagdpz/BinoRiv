% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function binoriv_percep_switch(w,superblock,triad)

global VAR
global SETTINGS

if VAR.eye_track == 1
    durTrialMax = VAR.trial_len; % max trial duration [s]
    fsMatlab = 1000; % matlab "sampling" frequency [Hz]
    bufferSize = durTrialMax * fsMatlab;
    memoryBuffer_fixloc = []; % [fix_loc_label_l fix_loc_label_r]
    memoryBuffer_eyepos = []; % [time trial x_eye y_eye]
    memoryBuffer_eyesize = []; % [time trial width hight]
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

% make binocular fixation spots so that they cannot overlap at one location
%
fix_loc_label = randi([1 4],2,VAR.num_trial);
if VAR.grating_task == 1
    for i=1:VAR.num_trial
        if fix_loc_label(1,i) == 1; while fix_loc_label(2,i) == 1 || fix_loc_label(2,i) == 4; fix_loc_label(2,i) = randi([1 4],1,1); end; end
        if fix_loc_label(1,i) == 2; while fix_loc_label(2,i) == 2 || fix_loc_label(2,i) == 3; fix_loc_label(2,i) = randi([1 4],1,1); end; end
        if fix_loc_label(1,i) == 3; while fix_loc_label(2,i) == 3 || fix_loc_label(2,i) == 1; fix_loc_label(2,i) = randi([1 4],1,1); end; end
        if fix_loc_label(1,i) == 4; while fix_loc_label(2,i) == 4 || fix_loc_label(2,i) == 2; fix_loc_label(2,i) = randi([1 4],1,1); end; end
    end
else
    for i=1:VAR.num_trial
        if fix_loc_label(1,i) == fix_loc_label(2,i)
            while fix_loc_label(1,i) == fix_loc_label(2,i)
                fix_loc_label(1,i) = randi([1 4],1,1);
                fix_loc_label(2,i) = randi([1 4],1,1);
            end
        end
    end
end


% fix_loc_label = randi([1 2],2,VAR.num_trial);
% if VAR.grating_task == 1
%     for i=1:VAR.num_trial
%         if fix_loc_label(1,i) == 1; fix_loc_label(2,i) = 2; end
%         if fix_loc_label(1,i) == 2; fix_loc_label(2,i) = 1; end
%     end
% else
%     for i=1:VAR.num_trial
%         if fix_loc_label(1,i) == fix_loc_label(2,i)
%             while fix_loc_label(1,i) == fix_loc_label(2,i)
%                 fix_loc_label(1,i) = randi([1 2],1,1);
%                 fix_loc_label(2,i) = randi([1 2],1,1);
%             end
%         end
%     end
% end

%% animation
[vbl, start] = Screen('Flip', w);
GetSecs;
% for online plotting
if VAR.eye_track == 1; figure('Name','Online GUI','Position',SETTINGS.GUI_coordinates,'Renderer', 'Painters'); axis([-5 5 -5 5]); axis square; grid on; hold on; end 
% main loop
for t = 1:VAR.num_trial
%     fix_loc_l = VAR.potential_loc(fix_loc_label(1,t),:);
    fix_loc_l = VAR.potential_loc(1,:);
    if VAR.grating_task == 1
%        fix_loc_r = VAR.potential_loc_grat(fix_loc_label(2,t),:);
        fix_loc_r = VAR.potential_loc_grat(4,:);
    else
        fix_loc_r = VAR.potential_loc(fix_loc_label(2,t),:);
    end
   
    reco = true;
    while ((GetSecs - start) < (t*VAR.trial_len))
        % fixation spots
        %Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        %Screen('FillOval', w, red, fix_loc_l);
        %Screen('FillOval', w, blue, fix_loc_r);
        %stimuli
        %Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [0 0 0 1]);
        %Screen('DrawTexture', w, imagetex_annulus); % annulus
        %Screen('DrawTexture', w, imagetex_bino);
        %Screen('DrawTexture', w, imagetex_l);
        %Screen('DrawTexture', w, imagetex_r);
        
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 1]);
        Screen('FrameOval', w, VAR.lineColour , VAR.ann_rect, VAR.lineWidth); % annulus
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 0 1 0]);
        Screen('DrawTexture', w, VAR.imagetex_bino);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 0 0 0]);
%         Screen('FillOval', w, VAR.red, fix_loc_l);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [0 0 1 0]);
%         Screen('FillOval', w, VAR.blue, fix_loc_r);
        
        % centre fixation spot 
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 0]);
        Screen('FillOval', w, [255 255 255], [VAR.rect(3)/2-VAR.radius VAR.rect(4)/2-VAR.radius VAR.rect(3)/2+VAR.radius VAR.rect(4)/2+VAR.radius]);

        %{
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 0 0 0]);
        Screen('DrawTexture', w, imagetex_l);
        Screen('FillOval', w, red, fix_loc_l);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [0 0 1 1]);
        Screen('DrawTexture', w, imagetex_r);
        Screen('FillOval', w, blue, fix_loc_r);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 1]);
        Screen('FrameOval', w, lineColour , ann_rect, lineWidth); % annulus
        %}
        vbl = Screen('Flip', w); % update screen
            
        [keyIsDown, press, KeyCode] = KbCheck;
        if KeyCode(VAR.EscapeKey)==1 
            break
        end
        
        while reco
            [keyIsDown, press, KeyCode] = KbCheck;
            if  KeyCode(VAR.EscapeKey)==1 
                break
            end
            reco_inner_ite =true;
            GetSecs;
            while reco_inner_ite
                time_sample = GetSecs - start; % sec
                if VAR.report == 1
                    % Record the time if the subject presses button
                    if KeyCode(VertKey)==1 && pre_VertKey==0
                        pre_VertKey = 1;
                        vert_press = vertcat(vert_press, time_sample*1000); % ms
                    elseif KeyCode(VertKey)==0 && pre_VertKey==1
                        pre_VertKey = 0;
                        vert_press = vertcat(vert_press, time_sample*1000); % ms
                    elseif KeyCode(HorzKey)==1 && pre_HorzKey==0
                        pre_HorzKey = 1;
                        horz_press = vertcat(horz_press, time_sample*1000); % ms
                    elseif KeyCode(HorzKey)==0 && pre_HorzKey==1
                        pre_HorzKey = 0;
                        horz_press = vertcat(horz_press, time_sample*1000); % ms
                    end
                end
                if VAR.eye_track == 1
                    [x_eye, y_eye] = aux_GetCalibratedEyeHandPos();
                    [width,hight] = vpx_GetPupilSize();
                    counter = counter + 1;
                    %memoryBuffer_fixloc(counter,1:4) = fix_loc_l;
                    %memoryBuffer_fixloc(counter,5:8) = fix_loc_r;
                    memoryBuffer_fixloc(counter,1) = fix_loc_label(1,t);
                    memoryBuffer_fixloc(counter,2) = fix_loc_label(2,t);
                    memoryBuffer_eyepos(counter,1) = time_sample*1000;
                    memoryBuffer_eyepos(counter,2) = t;
                    memoryBuffer_eyepos(counter,3) = x_eye; % in deg
                    memoryBuffer_eyepos(counter,4) = y_eye;
                    memoryBuffer_eyesize(counter,1) = time_sample*1000;
                    memoryBuffer_eyesize(counter,2) = t;
                    memoryBuffer_eyesize(counter,3) = width;
                    memoryBuffer_eyesize(counter,4) = hight;

                    % Online plot % heavy
    %                 % from px to deg|fixation spot
    %                 fix_loc_l_xcent = fix_loc_l(1,1)+(fix_loc_l(1,3)-fix_loc_l(1,1)) - VAR.rect(3)/2; % px from centre
    %                 fix_loc_l_ycent = -(fix_loc_l(1,2)+(fix_loc_l(1,4)-fix_loc_l(1,2)) - VAR.rect(4)/2);
    %                 fix_loc_r_xcent = fix_loc_r(1,1)+(fix_loc_r(1,3)-fix_loc_r(1,1)) - VAR.rect(3)/2; % px from centre
    %                 fix_loc_r_ycent = -(fix_loc_r(1,2)+(fix_loc_r(1,4)-fix_loc_r(1,2)) - VAR.rect(4)/2);
    %                 width_cmperpx = SETTINGS.screen_w_cm/VAR.rect(3); % cm/px of the screen
    %                 fix_loc_l_xcent = fix_loc_l_xcent*width_cmperpx; % cm from centre
    %                 fix_loc_l_ycent = fix_loc_l_ycent*width_cmperpx;
    %                 fix_loc_r_xcent = fix_loc_r_xcent*width_cmperpx; % cm from centre
    %                 fix_loc_r_ycent = fix_loc_r_ycent*width_cmperpx;
    %                 fix_loc_l_xcent = 2*atand(fix_loc_l_xcent/(2*SETTINGS.vd)); % deg
    %                 fix_loc_l_ycent = 2*atand(fix_loc_l_ycent/(2*SETTINGS.vd));
    %                 fix_loc_r_xcent = 2*atand(fix_loc_r_xcent/(2*SETTINGS.vd)); % deg
    %                 fix_loc_r_ycent = 2*atand(fix_loc_r_ycent/(2*SETTINGS.vd));
    %                 try 
    %                     x_eye_pre = memoryBuffer_eyepos(counter-1,3);
    %                     y_eye_pre = memoryBuffer_eyepos(counter-1,4);
    %                     line([x_eye_pre x_eye], [y_eye_pre y_eye], 'Color','k','Marker','none','LineWidth',0.2)
    %                     plot(x_eye,y_eye,'o', 'MarkerEdgeColor',[0.5 0.5 0.5])
    %                     plot(fix_loc_l_xcent,fix_loc_l_ycent,'*', 'MarkerEdgeColor','r')
    %                     plot(fix_loc_r_xcent,fix_loc_r_ycent,'*', 'MarkerEdgeColor','b')
    %                 catch
    %                 end
    %                 drawnow
    
        %         if VAR.eye_track == 1 % Very heavy, don't put in main loop
        %             filename = [VAR.fig_dir '/reco_online_bino_' num2str(superblock) '_' num2str(triad) '_' num2str(t) '.png'];
        %             saveas(gcf,filename)
        %             filename = [VAR.fig_dir '/reco_online_bino_' num2str(superblock) '_' num2str(triad) '_' num2str(t) '.fig'];
        %             saveas(gcf,filename)
        %         end
                end
                WaitSecs('Untiltime',GetSecs+SETTINGS.ite); % maintain matlab sampling rate
                reco_inner_ite = false;
            end
            if time_sample > t*VAR.trial_len; reco=false; end
        end
    end
end

%% save
if VAR.report == 1
    vert_press = array2table(vert_press);
    horz_press = array2table(horz_press);
    filename = [VAR.subj_dist '/report/bino/vertpress_repo_' num2str(superblock) '_' num2str(triad) '.csv'];
    writetable(vert_press, filename, 'WriteVariableNames', false);
    filename = [VAR.subj_dist '/report/bino/horzpress_repo_' num2str(superblock) '_' num2str(triad) '.csv'];
    writetable(horz_press, filename, 'WriteVariableNames', false);
end
    
if VAR.eye_track == 1
    if VAR.report == 1
        filename = [VAR.repo_dir '/bino/fixloc_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_fixloc);
        writetable(buffer, filename, 'WriteVariableNames', false);
        filename = [VAR.repo_dir '/bino/eyepos_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_eyepos);
        writetable(buffer, filename, 'WriteVariableNames', false);
        filename = [VAR.repo_dir '/bino/eyesize_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_eyesize);
        writetable(buffer, filename, 'WriteVariableNames', false);
    else
        filename = [VAR.norepo_dir '/bino/fixloc_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_fixloc);
        writetable(buffer, filename, 'WriteVariableNames', false);
        filename = [VAR.norepo_dir '/bino/eyepos_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_eyepos);
        writetable(buffer, filename, 'WriteVariableNames', false);
        filename = [VAR.norepo_dir '/bino/eyesize_' num2str(superblock) '_' num2str(triad) '.csv'];
        buffer = array2table(memoryBuffer_eyesize);
        writetable(buffer, filename, 'WriteVariableNames', false);
    end
end
