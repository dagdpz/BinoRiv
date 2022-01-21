% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function binoriv_percep_switch(w,red,blue,superblock,triad,trial_len,num_trial,imagetex_l,imagetex_r,imagetex_bino,potential_loc,report,subj_dist,EscapeKey,lineWidth, ann_rect, lineColour)

if report == 1
    % Define the key for reporting the content of perception
    KbName('UnifyKeyNames');
    EscapeKey = KbName('ESCAPE');
    VertKey = KbName('UpArrow');
    HorzKey = KbName('LeftArrow');
end

% make binocular fixation spots so that they cannot overlap at one location
fix_loc_label = randi([1 4],2,num_trial);
for i=1:num_trial
    if fix_loc_label(1,i) == fix_loc_label(2,i)
        while fix_loc_label(1,i) == fix_loc_label(2,i)
            fix_loc_label(1,i) = randi([1 4],1,1);
            fix_loc_label(2,i) = randi([1 4],1,1);
        end
    end
end

% animation
[vbl, start] = Screen('Flip', w);
for t = 1:num_trial
%    [vbl, start] = Screen('Flip', w); % do not Flip here; otherwise the screen will flicker.
    fix_loc_l = potential_loc(fix_loc_label(1,t),:);
    fix_loc_r = potential_loc(fix_loc_label(2,t),:); 
    if report == 1
        vert_press = [];
        horz_press = [];
    end
    
    while (vbl < (t*trial_len + start))
        % fixation spots
        %Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        %Screen('FillOval', w, red, fix_loc_l);
        %Screen('FillOval', w, blue, fix_loc_r);
        %stimuli
        %Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [0 0 0 1]);
        %Screen('DrawTexture', w, imagetex_annulus); % annulus
        %Screen('DrawTexture', w, imagetex_bino);
        %Screen('DrawTexture', w, imagetex_l);
%        Screen('DrawTexture', w, imagetex_r);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 0 0 0]);
        Screen('DrawTexture', w, imagetex_l);
        Screen('FillOval', w, red, fix_loc_l);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [0 0 1 1]);
        Screen('DrawTexture', w, imagetex_r);
        Screen('FillOval', w, blue, fix_loc_r);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 1]);
        Screen('FrameOval', w, lineColour , ann_rect, lineWidth); % annulus
             
        vbl = Screen('Flip', w); % update screen
        
        [keyIsDown, press, KeyCode] = KbCheck;
        if KeyCode(EscapeKey)==1 
            break
        end
        if report == 1
            % Record the time if subject pressed button
             if KeyCode(VertKey)==1
                vert_press = vertcat(vert_press, (vbl-start)*1000); % ms
            elseif KeyCode(HorzKey)==1 
                horz_press = vertcat(horz_press, (vbl-start)*1000); % ms
            end
        end
    end
end
if report == 1
    vert_press = array2table(vert_press);
    horz_press = array2table(horz_press);
    filename = [subj_dist '/report/bino/vertpress_repo_' num2str(superblock) '_' num2str(triad) '.csv'];
    writetable(vert_press, filename, 'WriteVariableNames', false);
    filename = [subj_dist '/report/bino/horzpress_repo_' num2str(superblock) '_' num2str(triad) '.csv'];
    writetable(horz_press, filename, 'WriteVariableNames', false);
end
