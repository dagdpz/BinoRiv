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

%% animation
[vbl, start] = Screen('Flip', w);
for t = 1:VAR.num_trial
    fix_loc_l = VAR.potential_loc(fix_loc_label(1,t),:);
    if VAR.grating_task == 1
        fix_loc_r = VAR.potential_loc_grat(fix_loc_label(2,t),:);
    else
        fix_loc_r = VAR.potential_loc(fix_loc_label(2,t),:);
    end
   
    while (vbl < (t*VAR.trial_len + start))
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
        
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 0 1 0]);
        Screen('DrawTexture', w, VAR.imagetex_bino);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 0 0 0]);
        Screen('FillOval', w, VAR.red, fix_loc_l);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [0 0 1 0]);
        Screen('FillOval', w, VAR.blue, fix_loc_r);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 1]);
        Screen('FrameOval', w, VAR.lineColour , VAR.ann_rect, VAR.lineWidth); % annulus

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
            %memoryBuffer_fixloc(counter,1:4) = fix_loc_l;
            %memoryBuffer_fixloc(counter,5:8) = fix_loc_r;
            memoryBuffer_fixloc(counter,1) = fix_loc_label(1,t);
            memoryBuffer_fixloc(counter,2) = fix_loc_label(2,t);
            memoryBuffer_eyepos(counter,1) = (vbl-start)*1000;
            memoryBuffer_eyepos(counter,2) = t;
            memoryBuffer_eyepos(counter,3) = x_eye; % in deg
            memoryBuffer_eyepos(counter,4) = y_eye;
            memoryBuffer_eyesize(counter,1) = (vbl-start)*1000;
            memoryBuffer_eyesize(counter,2) = t;
            memoryBuffer_eyesize(counter,3) = width;
            memoryBuffer_eyesize(counter,4) = hight;

            %{
            % Open graphics ("GUIs")
            GUI_fig_handle = figure('Name','Online GUI','Position',SETTINGS.GUI_coordinates,'Renderer', 'Painters');
            set(gca,'Ylim',[-SETTINGS.screen_lh_deg SETTINGS.screen_uh_deg],'Xlim',[-SETTINGS.screen_w_deg/2 SETTINGS.screen_w_deg/2],'XLimMode','manual','YLimMode','manual');
            if SETTINGS.matlab_version(1)>=2014
                dyn.eye_position_handle = line(-1000,-1000,'Color','r','Marker','o','XlimInclude','off','YLimInclude','off');%'background');
                dyn.hnd_position_handle = line(-1000,-1000,'Color','g','Marker','o','XlimInclude','off','YLimInclude','off');%,'background');
            else
                dyn.eye_position_handle = line(-1000,-1000,'Color','r','Marker','o','XlimInclude','off','YLimInclude','off','EraseMode','xor');%'background');
                dyn.hnd_position_handle = line(-1000,-1000,'Color','g','Marker','o','XlimInclude','off','YLimInclude','off','EraseMode','xor');%,'background');
            end
            hold on; drawnow;

            % Draw GUI
            dyn.correct_choice_target = trial(dyn.trialNumber).task.correct_choice_target;
            dyn.hnd_color=nanmean([task.hnd.fix.color_bright],1);
            if SETTINGS.GUI
                set(0, 'currentfigure', GUI_fig_handle);
                ang=0:0.02:2*pi;
                %aux_draw_GUI_targets(dyn,trial,eff,pha,ang)
                aux_draw_GUI_targets(dyn,trial,'eye','fix',ang)
                aux_draw_GUI_targets(dyn,trial,'eye','fi2',ang)
                aux_draw_GUI_targets(dyn,trial,'eye','tar',ang)
                aux_draw_GUI_targets(dyn,trial,'eye','ta2',ang)
                aux_draw_GUI_targets(dyn,trial,'eye','cue',ang)
                aux_draw_GUI_targets(dyn,trial,'hnd','fix',ang)
                aux_draw_GUI_targets(dyn,trial,'hnd','fi2',ang)
                aux_draw_GUI_targets(dyn,trial,'hnd','tar',ang)
                aux_draw_GUI_targets(dyn,trial,'hnd','ta2',ang)
                aux_draw_GUI_targets(dyn,trial,'hnd','cue',ang)
                if any(~isnan(task.hnd.fix.color_bright))
                    set(dyn.hnd_position_handle,'Color',nanmean([task.hnd.fix.color_bright],1)/255);
                end
            end
            %}
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
    filename = [VAR.subj_dist '/report/bino/fixloc_' num2str(superblock) '_' num2str(triad) '.csv'];
    buffer = array2table(memoryBuffer_fixloc);
    writetable(buffer, filename, 'WriteVariableNames', false);
    filename = [VAR.subj_dist '/report/bino/eyepos_' num2str(superblock) '_' num2str(triad) '.csv'];
    buffer = array2table(memoryBuffer_eyepos);
    writetable(buffer, filename, 'WriteVariableNames', false);
    filename = [VAR.subj_dist '/report/bino/eyesize_' num2str(superblock) '_' num2str(triad) '.csv'];
    buffer = array2table(memoryBuffer_eyesize);
    writetable(buffer, filename, 'WriteVariableNames', false);
end
