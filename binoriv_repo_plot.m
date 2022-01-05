% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function binoriv_repo_plot(subj,num_superblock,num_triad,num_trial,trial_len)

close all
subj = 'Ryo'
num_superblock = 1
num_triad = 1
num_trial = 8
trial_len = 2

subj_dir = fullfile('recording/', subj);
repo_dir = [subj_dir '/report'];
phys_dir = [subj_dir '/report/phys'];
bino_dir = [subj_dir '/report/bino'];

mkdir(repo_dir, 'figures')
fig_dir = [repo_dir '/figures'];

f1 = figure;
f2 = figure;

% Physical
figure(f1)
for spr = 1:num_superblock
    for trd = 1:num_triad
        vert_triad = [];
        horz_triad = [];
        ansfile = readtable([phys_dir '/answer' num2str(spr) '_' num2str(trd) '.csv'], 'ReadVariableNames', false);
        for trl = 1:num_trial
            vertfile = readtable([phys_dir '/vertpress_repo_' num2str(spr) '_' num2str(trd) '_' num2str(trl) '.csv']);
            %for i = 1:height(vertfile); xline(vertfile{i,1} + trial_len*1000*(trl-1), 'color', [0 0.4470 0.7410]); end % Blue
            for i = 1:height(vertfile); vert_triad = vertcat(vert_triad, vertfile{i,1} + trial_len*1000*(trl-1)); end
            horzfile = readtable([phys_dir '/horzpress_repo_' num2str(spr) '_' num2str(trd) '_' num2str(trl) '.csv']);
            %for i = 1:height(horzfile); xline(horzfile{i,1} + trial_len*1000*(trl-1), 'color', [0.8500 0.3250 0.0980]); end % Orange
            for i = 1:height(horzfile); horz_triad = vertcat(horz_triad, horzfile{i,1} + trial_len*1000*(trl-1)); end
        end
        % main plot
        try
            for i = 1:height(vert_triad)
                if (vert_triad(i+1,1)-vert_triad(i,1)) < 50
                    Square_colouring([vert_triad(i,1) vert_triad(i+1,1)], [0 0.4470 0.7410], [1 1], 0);
                end
            end
        catch
        end
        try
            for i = 1:height(horz_triad)
                if (horz_triad(i+1,1)-horz_triad(i,1)) < 50
                    Square_colouring([horz_triad(i,1) horz_triad(i+1,1)], [0.8500 0.3250 0.0980], [1 1], 0);
                end
            end
        catch
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
        xlabel('Response time [ms]')
        xlim([0 1000*trial_len*trl])
        yticks([]) % delete y-axis tick values
        %legend('Vertical', 'Horizontal', '', '', '', '', '', '', '', '')

        filename = [fig_dir '/phys_verthorzPress_' num2str(spr) '_' num2str(trd) '.png'];
        saveas(gcf,filename)

        hold on
        
    end
end

% BinoRiv
figure(f2)
for spr = 1:num_superblock
    for trd = 1:num_triad
        vert_triad = [];
        horz_triad = [];
        for trl = 1:num_trial
            vertfile = readtable([bino_dir '/vertpress_repo_' num2str(spr) '_' num2str(trd) '_' num2str(trl) '.csv']);
            %for i = 1:height(vertfile); xline(vertfile{i,1} + trial_len*1000*(trl-1), 'color', [0 0.4470 0.7410]); end % Blue
            for i = 1:height(vertfile); vert_triad = vertcat(vert_triad, vertfile{i,1} + trial_len*1000*(trl-1)); end
            horzfile = readtable([bino_dir '/horzpress_repo_' num2str(spr) '_' num2str(trd) '_' num2str(trl) '.csv']);
            %for i = 1:height(horzfile); xline(horzfile{i,1} + trial_len*1000*(trl-1), 'color', [0.8500 0.3250 0.0980]); end % Orange
            for i = 1:height(horzfile); horz_triad = vertcat(horz_triad, horzfile{i,1} + trial_len*1000*(trl-1)); end
        end
        % main plot
        try
            for i = 1:height(vert_triad)
                if (vert_triad(i+1,1)-vert_triad(i,1)) < 50
                    Square_colouring([vert_triad(i,1) vert_triad(i+1,1)], [0 0.4470 0.7410], [1 1], 0);
                end
            end
        catch
        end
        try
            for i = 1:height(horz_triad)
                if (horz_triad(i+1,1)-horz_triad(i,1)) < 50
                    Square_colouring([horz_triad(i,1) horz_triad(i+1,1)], [0.8500 0.3250 0.0980], [1 1], 0);
                end
            end
        catch
        end
        % draw border lines
        for bl = (trial_len*1000):(trial_len*1000):(trial_len*1000*(num_trial-1)); xline(bl, '--', 'LineWidth', 1.0); end

        xlabel('Response time [ms]')
        xlim([0 1000*trial_len*trl])
        yticks([]) % delete y-axis tick values
        %legend('Vertical', 'Horizontal', '', '', '', '', '', '', '', '')

        filename = [fig_dir '/bino_verthorzPress_' num2str(spr) '_' num2str(trd) '.png'];
        saveas(gcf,filename)

        hold on
        
    end
end