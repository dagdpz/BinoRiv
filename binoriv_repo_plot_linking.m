function binoriv_repo_plot_linking(subj,num_superblock,num_triad,num_trial,trial_len)

close all
subj_type = 0;
subj = 'Ryo'
num_superblock = 1
num_triad = 1
num_trial = 8
trial_len = 2
report = 1
colour_comb = 1; % 0 is (left:Red right:Blue), 1 is (left:Blue right:Red)
num_plot = 2

if subj_type == 0
    subj_dir = fullfile('recording/human/', subj);
elseif subj_type == 1
    subj_dir = fullfile('recording/monkey/', subj);
end
if report == 0
    filename = [subj_dir '/report/switch.mat'];
elseif report == 1
    filename = [subj_dir '/noreport/switch.mat'];
end
repo_dir = [subj_dir '/report'];
phys_dir = [subj_dir '/report/phys'];
bino_dir = [subj_dir '/report/bino'];

fig_dir = [repo_dir '/figures'];


for spb = 1:num_superblock
    % load switch log
    if report == 0
        switchfile = readtable([subj_dir '/noreport/switch_' num2str(spb) '.csv'], 'ReadVariableNames', false);
    elseif report == 1
        switchfile = readtable([subj_dir '/report/switch_' num2str(spb) '.csv'], 'ReadVariableNames', false);
    end
        
    for trd = 1:num_triad
        % load data
        ansfile = readtable([phys_dir '/answer_' num2str(spb) '_' num2str(trd) '.csv'], 'ReadVariableNames', false);
        
        vert_triad_phys = [];
        horz_triad_phys = [];
        vertfile_phys = readtable([phys_dir '/vertpress_repo_' num2str(spb) '_' num2str(trd) '.csv']);
        for i = 1:height(vertfile_phys); vert_triad_phys = vertcat(vert_triad_phys, vertfile_phys{i,1}); end
        horzfile_phys = readtable([phys_dir '/horzpress_repo_' num2str(spb) '_' num2str(trd) '.csv']);
        for i = 1:height(horzfile_phys); horz_triad_phys = vertcat(horz_triad_phys, horzfile_phys{i,1}); end
        vert_triad_br = [];
        horz_triad_br = [];
        vertfile_br = readtable([bino_dir '/vertpress_repo_' num2str(spb) '_' num2str(trd) '.csv']);
        for i = 1:height(vertfile_br); vert_triad_br = vertcat(vert_triad_br, vertfile_br{i,1}); end
        horzfile_br = readtable([bino_dir '/horzpress_repo_' num2str(spb) '_' num2str(trd) '.csv']);
        for i = 1:height(horzfile_br); horz_triad_br = vertcat(horz_triad_br, horzfile_br{i,1}); end
        
        % main plot
        %% phys->bino
        if switchfile{trd,1} == 0 % if 0 phys->bino, if 1 bino->phys
            subplot(num_plot,1,1);
            for trl = 1:num_trial
                answer = ansfile{trl,1};
                answer = cell2mat(answer);
                if answer == 'L'
                    plot([1000*trial_len*(trl-1) 1000*trial_len*(trl)], [2 2], 'color', 'k')
                    hold on
                elseif answer == 'R'
                    plot([1000*trial_len*(trl-1) 1000*trial_len*(trl)], [1 1], 'color', 'k')
                end
            end
            plot([1000*trial_len*num_trial 1000*trial_len*num_trial*2], [0 0], 'color', 'k') % BR
            xticks(0:1000*trial_len:1000*trial_len*num_trial*2)
            xlim([0 1000*trial_len*num_trial*2])
            ylim([-0.5 2.5])
            yticks([0 1 2])
            yticklabels({'BR','Blue','Red'})
            grid on
            
            subplot(num_plot,1,2);
            try % phys
                for i = 1:height(vert_triad_phys)
                    if rem(i,2) == 1; Square_colouring([vert_triad_phys(i,1) vert_triad_phys(i+1,1)], [0.8 0.2 0], [1 1], 0); end
                end
            catch
                Square_colouring([vert_triad_phys(i,1) 1000*trial_len*num_trial], [0.8 0.2 0], [1 1], 0);
            end
            try
                for i = 1:height(horz_triad_phys)
                    if rem(i,2) == 1; Square_colouring([horz_triad_phys(i,1) horz_triad_phys(i+1,1)], [0 0.2 0.8], [1 1], 0); end
                end
            catch
                Square_colouring([horz_triad_phys(i,1) 1000*trial_len*num_trial], [0 0.2 0.8], [1 1], 0);
            end
            try % BR
            for i = 1:height(vert_triad_br)
                if rem(i,2) == 1; Square_colouring([(1000*trial_len*num_trial)+vert_triad_br(i,1) (1000*trial_len*num_trial)+vert_triad_br(i+1,1)], [0.8 0.2 0], [1 1], 0); end
            end
            catch
                Square_colouring([(1000*trial_len*num_trial)+vert_triad_br(i,1) (1000*trial_len*num_trial)+1000*trial_len*num_trial], [0.8 0.2 0], [1 1], 0);
            end
            try
                for i = 1:height(horz_triad_br)
                    if rem(i,2) == 1; Square_colouring([(1000*trial_len*num_trial)+horz_triad_br(i,1) (1000*trial_len*num_trial)+horz_triad_br(i+1,1)], [0 0.2 0.8], [1 1], 0); end
                end
            catch
                Square_colouring([(1000*trial_len*num_trial)+horz_triad_br(i,1) (1000*trial_len*num_trial)+1000*trial_len*num_trial], [0 0.2 0.8], [1 1], 0);
            end
            % draw border lines
            for bl = (trial_len*1000):(trial_len*1000):(trial_len*1000*(num_trial*2-1)); xline(bl, '--', 'LineWidth', 1.0); end
            xlabel('Time [ms]')
            xticks(0:1000*trial_len:1000*trial_len*num_trial*2)
            xlim([0 1000*trial_len*num_trial*2])
            yticks([]) % delete y-axis tick values
            
            filename = [fig_dir '/reco_' num2str(spb) '_' num2str(trd) '.png'];
            saveas(gcf,filename)
            filename = [fig_dir '/reco_' num2str(spb) '_' num2str(trd) '.fig'];
            saveas(gcf,filename)
            clf

       %% bino->phys
        elseif switchfile{trd,1} == 1 % if 0 phys->bino, if 1 bino->phys
            subplot(num_plot,1,1);
            plot([0 1000*trial_len*num_trial], [0 0], 'color', 'k') % BR
            hold on
            for trl = 1:num_trial
                answer = ansfile{trl,1};
                answer = cell2mat(answer);
                if answer == 'L'
                    plot([(1000*trial_len*num_trial)+(1000*trial_len*(trl-1)) (1000*trial_len*num_trial)+(1000*trial_len*trl)], [2 2], 'color', 'k')
                elseif answer == 'R'
                    plot([(1000*trial_len*num_trial)+(1000*trial_len*(trl-1)) (1000*trial_len*num_trial)+(1000*trial_len*trl)], [1 1], 'color', 'k')
                end
            end
            xticks(0:1000*trial_len:1000*trial_len*num_trial*2)
            xlim([0 1000*trial_len*num_trial*2])
            ylim([-0.5 2.5])
            yticks([0 1 2])
            yticklabels({'BR','Blue','Red'})
            grid on
            
            subplot(num_plot,1,2);
            try % BR
            for i = 1:height(vert_triad_br)
                if rem(i,2) == 1; Square_colouring([vert_triad_br(i,1) vert_triad_br(i+1,1)], [0.8 0.2 0], [1 1], 0); end
            end
            catch
                Square_colouring([vert_triad_br(i,1) 1000*trial_len*num_trial], [0.8 0.2 0], [1 1], 0);
            end
            try
                for i = 1:height(horz_triad_br)
                    if rem(i,2) == 1; Square_colouring([horz_triad_br(i,1) horz_triad_br(i+1,1)], [0 0.2 0.8], [1 1], 0); end
                end
            catch
                Square_colouring([horz_triad_br(i,1) 1000*trial_len*num_trial], [0 0.2 0.8], [1 1], 0);
            end
            try % phys
                for i = 1:height(vert_triad_phys)
                    if rem(i,2) == 1; Square_colouring([(1000*trial_len*num_trial)+vert_triad_phys(i,1) (1000*trial_len*num_trial)+vert_triad_phys(i+1,1)], [0.8 0.2 0], [1 1], 0); end
                end
            catch
                Square_colouring([(1000*trial_len*num_trial)+vert_triad_phys(i,1) (1000*trial_len*num_trial)+1000*trial_len*num_trial], [0.8 0.2 0], [1 1], 0);
            end
            try
                for i = 1:height(horz_triad_phys)
                    if rem(i,2) == 1; Square_colouring([(1000*trial_len*num_trial)+horz_triad_phys(i,1) (1000*trial_len*num_trial)+horz_triad_phys(i+1,1)], [0 0.2 0.8], [1 1], 0); end
                end
            catch
                Square_colouring([(1000*trial_len*num_trial)+horz_triad_phys(i,1) (1000*trial_len*num_trial)+1000*trial_len*num_trial], [0 0.2 0.8], [1 1], 0);
            end
            % draw border lines
            for bl = (trial_len*1000):(trial_len*1000):(trial_len*1000*(num_trial*2-1)); xline(bl, '--', 'LineWidth', 1.0); end
            xlabel('Time [ms]')
            xticks(0:1000*trial_len:1000*trial_len*num_trial*2)
            xlim([0 1000*trial_len*num_trial*2])
            yticks([]) % delete y-axis tick values
            
            filename = [fig_dir '/reco_' num2str(spb) '_' num2str(trd) '.png'];
            saveas(gcf,filename)
            filename = [fig_dir '/reco_' num2str(spb) '_' num2str(trd) '.fig'];
            saveas(gcf,filename)
            clf
        end
    end
end                
        
%xlim([0 1000*trial_len*num_trial*2])
