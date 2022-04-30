% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function binoriv_repo_plot(subj,num_superblock,num_triad,num_trial,trial_len)

close all
subj_type = 0;
subj = 'Ryo'
num_superblock = 1
num_triad = 11
num_trial = 8
trial_len = 2
eye_track = 0; % 0: eye tracking off, 1: eye tracking on

if subj_type == 0
    subj_dir = fullfile('recording/human/', subj);
elseif subj_type == 1
    subj_dir = fullfile('recording/monkey/', subj);
end
repo_dir = [subj_dir '/report'];
phys_dir = [subj_dir '/report/phys'];
bino_dir = [subj_dir '/report/bino'];

mkdir(repo_dir, 'figures')
fig_dir = [repo_dir '/figures'];

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;

% BinoRiv eye-tracking
if eye_track == 1
    figure(f1)
    time = [];
    eyepos_x = [];
    eyepos_y = [];
    for spr = 1:num_superblock
        for trd = 1:num_triad
            eyeposfile = readtable([bino_dir '/eyepos_' num2str(spr) '_' num2str(trd) '.csv']);

            for i = 1:height(eyeposfile); time = vertcat(time, eyeposfile{i,1}); end
            for i = 1:height(eyeposfile); eyepos_x = vertcat(eyepos_x, eyeposfile{i,3}); end
            for i = 1:height(eyeposfile); eyepos_y = vertcat(eyepos_y, eyeposfile{i,4}); end
            plot(time, eyepos_x)
            hold on
            plot(time, eyepos_y)
            xlabel('Time [ms]')
            xlim([0 1000*trial_len*num_trial])
            ylim([-1 1])
            legend('x-coordinate of gaze', 'y-coordinate of gaze')

            filename = [fig_dir '/bino_eyepos_' num2str(spr) '_' num2str(trd) '.png'];
            saveas(gcf,filename)
        end
    end
end

% Physical eye-tracking
if eye_track == 1
    figure(f2)
    time = [];
    eyepos_x = [];
    eyepos_y = [];
    for spr = 1:num_superblock
        for trd = 1:num_triad
            eyeposfile = readtable([phys_dir '/eyepos_' num2str(spr) '_' num2str(trd) '.csv']);

            for i = 1:height(eyeposfile); time = vertcat(time, eyeposfile{i,1}); end
            for i = 1:height(eyeposfile); eyepos_x = vertcat(eyepos_x, eyeposfile{i,3}); end
            for i = 1:height(eyeposfile); eyepos_y = vertcat(eyepos_y, eyeposfile{i,4}); end
            plot(time, eyepos_x)
            hold on
            plot(time, eyepos_y)
            xlabel('Time [ms]')
            xlim([0 1000*trial_len*num_trial])
            ylim([-1 1])
            legend('x-coordinate of gaze', 'y-coordinate of gaze')

            filename = [fig_dir '/phys_eyepos_' num2str(spr) '_' num2str(trd) '.png'];
            saveas(gcf,filename)
        end
    end
end

% BinoRiv
if eye_track == 1
    figure(f3)
else
    figure(f1)
end
for spr = 1:num_superblock
    for trd = 1:num_triad
        vert_triad = [];
        horz_triad = [];
        vertfile = readtable([bino_dir '/vertpress_repo_' num2str(spr) '_' num2str(trd) '.csv']);
        %for i = 1:height(vertfile); xline(vertfile{i,1} + trial_len*1000*(trl-1), 'color', [0 0.4470 0.7410]); end % Blue
        for i = 1:height(vertfile); vert_triad = vertcat(vert_triad, vertfile{i,1}); end
        horzfile = readtable([bino_dir '/horzpress_repo_' num2str(spr) '_' num2str(trd) '.csv']);
        %for i = 1:height(horzfile); xline(horzfile{i,1} + trial_len*1000*(trl-1), 'color', [0.8500 0.3250 0.0980]); end % Orange
        for i = 1:height(horzfile); horz_triad = vertcat(horz_triad, horzfile{i,1}); end
        
        % main plot
        try
            for i = 1:height(vert_triad)
%                if rem(i,2) == 1; Square_colouring([vert_triad(i,1) vert_triad(i+1,1)], [0 0.4470 0.7410], [1 1], 0); end
                if rem(i,2) == 1; Square_colouring([vert_triad(i,1) vert_triad(i+1,1)], [0.8 0.2 0], [1 1], 0); end
            end
        catch
            %Square_colouring([vert_triad(i,1) 1000*trial_len*num_trial], [0 0.4470 0.7410], [1 1], 0);
            Square_colouring([vert_triad(i,1) 1000*trial_len*num_trial], [0.8 0.2 0], [1 1], 0);
        end
        try
            for i = 1:height(horz_triad)
%                if rem(i,2) == 1; Square_colouring([horz_triad(i,1) horz_triad(i+1,1)], [0.8500 0.3250 0.0980], [1 1], 0); end
                if rem(i,2) == 1; Square_colouring([horz_triad(i,1) horz_triad(i+1,1)], [0 0.2 0.8], [1 1], 0); end
            end
        catch
            %Square_colouring([horz_triad(i,1) 1000*trial_len*num_trial], [0.8500 2 0.0980], [1 1], 0);
            Square_colouring([horz_triad(i,1) 1000*trial_len*num_trial], [0 0.2 0.8], [1 1], 0);
        end
        % draw border lines
        for bl = (trial_len*1000):(trial_len*1000):(trial_len*1000*(num_trial-1)); xline(bl, '--', 'LineWidth', 1.0); end

        xlabel('Time [ms]')
        xlim([0 1000*trial_len*num_trial])
        yticks([]) % delete y-axis tick values
        %legend('Vertical', 'Horizontal', '')

        filename = [fig_dir '/bino_verthorzPress_' num2str(spr) '_' num2str(trd) '.png'];
        saveas(gcf,filename)
        filename = [fig_dir '/bino_verthorzPress_' num2str(spr) '_' num2str(trd) '.fig'];
        saveas(gcf,filename)
        
        clf
    end
end

% Physical
if eye_track == 1
    figure(f4)
else
    figure(f2)
end
for spr = 1:num_superblock
    for trd = 1:num_triad
        vert_triad = [];
        horz_triad = [];
        ansfile = readtable([phys_dir '/answer_' num2str(spr) '_' num2str(trd) '.csv'], 'ReadVariableNames', false);
        vertfile = readtable([phys_dir '/vertpress_repo_' num2str(spr) '_' num2str(trd) '.csv']);
        %for i = 1:height(vertfile); xline(vertfile{i,1} + trial_len*1000*(trl-1), 'color', [0 0.4470 0.7410]); end % Blue
        for i = 1:height(vertfile); vert_triad = vertcat(vert_triad, vertfile{i,1}); end
        horzfile = readtable([phys_dir '/horzpress_repo_' num2str(spr) '_' num2str(trd) '.csv']);
        %for i = 1:height(horzfile); xline(horzfile{i,1} + trial_len*1000*(trl-1), 'color', [0.8500 0.3250 0.0980]); end % Orange
        for i = 1:height(horzfile); horz_triad = vertcat(horz_triad, horzfile{i,1}); end
        
        % main plot
        try
            for i = 1:height(vert_triad)
%                if rem(i,2) == 1; Square_colouring([vert_triad(i,1) vert_triad(i+1,1)], [0 0.4470 0.7410], [1 1], 0); end
                if rem(i,2) == 1; Square_colouring([vert_triad(i,1) vert_triad(i+1,1)], [0.8 0.2 0], [1 1], 0); end
            end
        catch
            %Square_colouring([vert_triad(i,1) 1000*trial_len*num_trial], [0 0.4470 0.7410], [1 1], 0);
            Square_colouring([vert_triad(i,1) 1000*trial_len*num_trial], [0.8 0.2 0], [1 1], 0);
        end
        try
            for i = 1:height(horz_triad)
%                if rem(i,2) == 1; Square_colouring([horz_triad(i,1) horz_triad(i+1,1)], [0.8500 0.3250 0.0980], [1 1], 0); end
                if rem(i,2) == 1; Square_colouring([horz_triad(i,1) horz_triad(i+1,1)], [0 0.2 0.8], [1 1], 0); end
            end
        catch
            %Square_colouring([horz_triad(i,1) 1000*trial_len*num_trial], [0.8500 0.3250 0.0980], [1 1], 0);
            Square_colouring([horz_triad(i,1) 1000*trial_len*num_trial], [0 0.2 0.8], [1 1], 0);
        end
        % draw border lines
        for bl = (trial_len*1000):(trial_len*1000):(trial_len*1000*(num_trial-1)); xline(bl, '--', 'LineWidth', 1.0); end
        % write answer labels
        for i = 1:height(ansfile)
            if cell2mat(ansfile{i,1}) == 'L'
                text(trial_len*1000*(i-1), 0.05, 'Vertical')
            elseif cell2mat(ansfile{i,1}) == 'R'
                text(trial_len*1000*(i-1), 0.05, 'Horizontal')
            end
        end
        xlabel('Time [ms]')
        xlim([0 1000*trial_len*num_trial])
        yticks([]) % delete y-axis tick values
        %legend('Vertical', 'Horizontal', '', '', '', '', '', '', '', '')

        filename = [fig_dir '/phys_verthorzPress_' num2str(spr) '_' num2str(trd) '.png'];
        saveas(gcf,filename)
        filename = [fig_dir '/phys_verthorzPress_' num2str(spr) '_' num2str(trd) '.fig'];
        saveas(gcf,filename)
        clf
    end
end


%{
% triad
figure(f5)
subplot(num_plot,1,1);
try
    for i = 1:height(vert_triad)
        if rem(i,2) == 1; Square_colouring([vert_triad(i,1) vert_triad(i+1,1)], [0 0.4470 0.7410], [1 1], 0); end
    end
catch
    Square_colouring([vert_triad(i,1) 1000*trial_len*num_trial], [0 0.4470 0.7410], [1 1], 0);
end
try
    for i = 1:height(horz_triad)
        if rem(i,2) == 1; Square_colouring([horz_triad(i,1) horz_triad(i+1,1)], [0.8500 0.3250 0.0980], [1 1], 0); end
    end
catch
    Square_colouring([horz_triad(i,1) 1000*trial_len*num_trial], [0.8500 0.3250 0.0980], [1 1], 0);
end
xlabel('Time [ms]')
xlim([0 1000*trial_len*num_trial*2])
%}