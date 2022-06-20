% Output: Success rate of gaze fixation
% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function gaze_fix_calc()

close all

%% load data
subj_type = 0;
subj = 'Rebecca_dc_test_0617';
num_superblock = 2;
num_triad = 12;
num_trial = 8;
trial_len = 2000; % ms
report = 0;
only_grating = 1;
colour_comb = 1; % 0 is (left:Red right:Blue), 1 is (left:Blue right:Red)
monitor_x = 2560;
monitor_y = 1440;
screen_width = 61-1.4;
dist_scr = 47; % distance from screen (cm)
saccade_allow = 350; % ms

threshold = 1; % deg, range you allow as a gaze fixation

if subj_type == 0
    subj_dir = fullfile('recording/human/', subj);
elseif subj_type == 1
    subj_dir = fullfile('recording/monkey/', subj);
end
if report == 0
    repo_dir = [subj_dir '/noreport'];
    phys_dir = [subj_dir '/noreport/phys'];
    bino_dir = [subj_dir '/noreport/bino'];
else
    repo_dir = [subj_dir '/report'];
    phys_dir = [subj_dir '/report/phys'];
    bino_dir = [subj_dir '/report/bino'];
end

fig_dir = [repo_dir '/figures'];
mkdir(fig_dir)

% load fixation spot location
% variables = load([repo_dir '/variables_repo_*.mat']);
variables = load([repo_dir '/variables_repo_20220617']);
fp_loc = variables.VAR.potential_loc;
try fp_loc_grat = variables.VAR.potential_loc_grat; catch; end

%% main
success = 0;
failure = 0;
result = [];
for spb = 1:num_superblock
    for trd = 1:num_triad
        % phys
        time = [];
        eyepos_x = [];
        eyepos_y = [];
        fix_loc = [];
        distance = [];
        eyeposfile_phys = readtable([phys_dir '/eyepos_' num2str(spb) '_' num2str(trd) '.csv']);
        fixlocfile_phys = readtable([phys_dir '/fixloc_' num2str(spb) '_' num2str(trd) '.csv']);
        for i = 1:height(eyeposfile_phys); time = vertcat(time, eyeposfile_phys{i,1}); end
        for i = 1:height(eyeposfile_phys); eyepos_x = vertcat(eyepos_x, eyeposfile_phys{i,3}); end %deg
        for i = 1:height(eyeposfile_phys); eyepos_y = vertcat(eyepos_y, eyeposfile_phys{i,4}); end
        for i = 1:height(fixlocfile_phys); fix_loc(i,1) = fixlocfile_phys{i,1}; end
        for i = 1:height(eyeposfile_phys)
            if fix_loc(i,1) == 1
                fix_loc_x = fp_loc(1,1)+(fp_loc(1,3) - fp_loc(1,1))/2; % px
                fix_loc_y = fp_loc(1,2)+(fp_loc(1,4) - fp_loc(1,2))/2;
            elseif fix_loc(i,1) == 2
                fix_loc_x = fp_loc_grat(4,1)+(fp_loc_grat(4,3) - fp_loc_grat(4,1))/2;
                fix_loc_y = fp_loc_grat(4,2)+(fp_loc_grat(4,4) - fp_loc_grat(4,2))/2;
            end
%             if only_grating == 0
%                 if fix_loc(i,1) == 1 % fix_left/upleft
%                     fix_loc_x = fp_loc(1,1)+(fp_loc(1,3) - fp_loc(1,1))/2; % px
%                     fix_loc_y = fp_loc(1,2)+(fp_loc(1,4) - fp_loc(1,2))/2;
%                 elseif fix_loc(i,1) == 2 % fix_right/upright
%                     fix_loc_x = fp_loc(2,1)+(fp_loc(2,3) - fp_loc(2,1))/2;
%                     fix_loc_y = fp_loc(2,2)+(fp_loc(2,4) - fp_loc(2,2))/2;
%                 elseif fix_loc(i,1) == 3 % fix_up/belowleft
%                     fix_loc_x = fp_loc(3,1)+(fp_loc(3,3) - fp_loc(3,1))/2;
%                     fix_loc_y = fp_loc(3,2)+(fp_loc(3,4) - fp_loc(3,2))/2;
%                 elseif fix_loc(i,1) == 4 % fix_below/belowright
%                     fix_loc_x = fp_loc(4,1)+(fp_loc(4,3) - fp_loc(4,1))/2;
%                     fix_loc_y = fp_loc(4,2)+(fp_loc(4,4) - fp_loc(4,2))/2;
%                 end
%             elseif only_grating == 1
%                 if fix_loc(i,1) == 1 % fix_left/upleft
%                     fix_loc_x = fp_loc_grat(1,1)+(fp_loc_grat(1,3) - fp_loc_grat(1,1))/2; % px
%                     fix_loc_y = fp_loc_grat(1,2)+(fp_loc_grat(1,4) - fp_loc_grat(1,2))/2;
%                 elseif fix_loc(i,1) == 2 % fix_right/upright
%                     fix_loc_x = fp_loc_grat(2,1)+(fp_loc_grat(2,3) - fp_loc_grat(2,1))/2;
%                     fix_loc_y = fp_loc_grat(2,2)+(fp_loc_grat(2,4) - fp_loc_grat(2,2))/2;
%                 elseif fix_loc(i,1) == 3 % fix_up/belowleft
%                     fix_loc_x = fp_loc_grat(3,1)+(fp_loc_grat(3,3) - fp_loc_grat(3,1))/2;
%                     fix_loc_y = fp_loc_grat(3,2)+(fp_loc_grat(3,4) - fp_loc_grat(3,2))/2;
%                 elseif fix_loc(i,1) == 4 % fix_below/belowright
%                     fix_loc_x = fp_loc_grat(4,1)+(fp_loc_grat(4,3) - fp_loc_grat(4,1))/2;
%                     fix_loc_y = fp_loc_grat(4,2)+(fp_loc_grat(4,4) - fp_loc_grat(4,2))/2;
%                 end
%             end
            
            % from px to deg|fixation spot
            fix_loc_x = fix_loc_x - monitor_x/2; % px from centre
            fix_loc_y = -(fix_loc_y - monitor_y/2);
            width_cmperpx = screen_width/monitor_x; % cm/px of the screen
            fix_loc_x = fix_loc_x * width_cmperpx; % cm from centre
            fix_loc_y = fix_loc_y * width_cmperpx;
            % V(deg) = 2*arctan(size/(2*DistFromMonitor))
            fix_loc_x = 2*atand(fix_loc_x/(2*dist_scr)); % deg
            fix_loc_y = 2*atand(fix_loc_y/(2*dist_scr));

            distance(i,1) = sqrt((eyepos_x(i,1)-fix_loc_x).^2+(eyepos_y(i,1)-fix_loc_y).^2);
            
%             figure('color','white');
%             % plot fixation spots
%             scatter(fix_loc_x, fix_loc_y, 100, 'red', '*'); hold on
%             % plot circle
%             t = linspace(0,2*pi,100);
%             cx = 0; cy = 0; % centre
%             r = 7.5/2;           % radius
%             plot(r*sin(t)+cx,r*cos(t)+cy)
%             axis square 
%             xlim([-4 4])
%             ylim([-4 4])
%             grid on
        end
       

        for trl = 1:num_trial            
            success_count = 0;
            if trl == 1
                num_saccade_sample = time < saccade_allow;
                num_saccade_sample = nnz(num_saccade_sample); % nnz counts non-zero values
                num_saccade_sample_post = time < (trl*trial_len)+saccade_allow;
                num_saccade_sample_post = nnz(num_saccade_sample_post);
            elseif trl == length(num_trial)
                num_saccade_sample = time < ((trl-1)*trial_len)+saccade_allow;
                num_saccade_sample = nnz(num_saccade_sample);
                num_saccade_sample_post = time < (trl*trial_len);
                num_saccade_sample_post = nnz(num_saccade_sample_post);
            else 
                num_saccade_sample = time < ((trl-1)*trial_len)+saccade_allow;
                num_saccade_sample = nnz(num_saccade_sample);
                num_saccade_sample_post = time < (trl*trial_len)+saccade_allow;
                num_saccade_sample_post = nnz(num_saccade_sample_post);
            end

            for i = num_saccade_sample:num_saccade_sample_post
                if distance(i,:) < threshold; success_count = success_count +1; end
            end
            
            if success_count < (num_saccade_sample_post-num_saccade_sample)/2
                sprintf('Fixation [superblock %d, triad %d, trial %d]: Failed', spb, trd, trl)
                failure = failure + 1;
                result = vertcat(result, "failed");
            elseif success_count > (num_saccade_sample_post-num_saccade_sample)/2
                sprintf('Fixation [superblock %d, triad %d, trial %d]: Successed', spb, trd, trl)
                success = success + 1;
                result = vertcat(result, "successed");
            end
            
            
        end

        % bino
    end
end
sprintf('Num of success: %d', success)
sprintf('Num of failure: %d', failure)