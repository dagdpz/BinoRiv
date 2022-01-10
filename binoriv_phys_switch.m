% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function phys_switch(w,red,blue,superblock,triad,trial_len,num_trial,imagetex_l,imagetex_r,potential_loc,report,phys_stim,subj_dist)

if report == 1
    % Define the key for reporting the content of perception
    KbName('UnifyKeyNames');
    EscapeKey = KbName('ESCAPE');
    VertKey = KbName('UpArrow');
    HorzKey = KbName('LeftArrow');
end
    
fix_loc_label = randi([1 4],1,num_trial);
answer = [];

% animation
for t = 1:num_trial
    [vbl, start] = Screen('Flip', w);
    fix_loc = potential_loc(fix_loc_label(1,t),:);
    if report == 1
        vert_press = [];
        horz_press = [];
        if phys_stim(t) == 0 % L
            answer = vertcat(answer, 'L');
        elseif phys_stim(t) == 1 % R
            answer = vertcat(answer, 'R');
        end
    end
    
    while (vbl < (trial_len + start))
        if phys_stim(t) == 0 % L
            Screen('DrawTexture', w, imagetex_l);
            Screen('FillOval', w, red, fix_loc)
        elseif phys_stim(t) == 1 % R
            Screen('DrawTexture', w, imagetex_r);
            Screen('FillOval', w, blue, fix_loc)
        end
        vbl = Screen('Flip', w); % return current time
        
        if report == 1
            % Record the time if subject pressed button
            [keyIsDown, press, KeyCode] = KbCheck;
            if KeyCode(EscapeKey)==1 
                break
            elseif KeyCode(VertKey)==1
                vert_press = vertcat(vert_press, (vbl-start)*1000); % ms
            elseif KeyCode(HorzKey)==1 
                horz_press = vertcat(horz_press, (vbl-start)*1000); % ms
            end
        end
    end
    
    if report == 1
        vert_press = array2table(vert_press);
        horz_press = array2table(horz_press);
        filename = [subj_dist '/report/phys/vertpress_repo_' num2str(superblock) '_' num2str(triad) '_' num2str(t) '.csv'];
        writetable(vert_press, filename, 'WriteVariableNames', false);
        filename = [subj_dist '/report/phys/horzpress_repo_' num2str(superblock) '_' num2str(triad) '_' num2str(t) '.csv'];
        writetable(horz_press, filename, 'WriteVariableNames', false);
    end
end

if report == 1
    answer = array2table(answer);
    filename = [subj_dist '/report/phys/answer' num2str(superblock) '_' num2str(triad) '.csv'];
    writetable(answer, filename, 'WriteVariableNames', false);
end