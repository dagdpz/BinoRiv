% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function binoriv_repo_plot_linking_ver2(subj,num_superblock,num_triad,num_trial,trial_len)

close all
subj_type = 0;
subj = 'Luba_0419'
num_superblock = 1
num_triad = 12
num_trial = 8
trial_len = 2
report = 1
only_grating = 1;
eye_track = 1;
colour_comb = 1; % 0 is (left:Red right:Blue), 1 is (left:Blue right:Red)
monitor_x = 2560;
monitor_y = 1440;
screen_width = 61-1.4;
dist_scr = 47; % distance from screen (cm)

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
mkdir(fig_dir)

ratestorage_br_vert = 0;
ratestorage_br_horz = 0;
ratestorage_br_mix = 0;
ratestorage_phys_vert = 0;
ratestorage_phys_horz = 0;

% load fixation spot location
%variables = load([repo_dir '/variables_repo_*.mat']);
variables = load([repo_dir '/variables_repo_20220419']);
fp_loc = variables.VAR.potential_loc;
try fp_loc_grat = variables.VAR.potential_loc_grat; catch; end


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
        
        timestorage_br_vert = 0;
        timestorage_br_horz = 0;
        timestorage_phys_vert = 0;
        timestorage_phys_horz = 0;             
        
        % main plot
        %% phys->bino
        if switchfile{trd,1} == 0 % if 0 phys->bino, if 1 bino->phys
            t = tiledlayout(1,1);
            ax1 = axes(t);
            if eye_track == 1; yyaxis left; end
            for trl = 1:num_trial
                answer = ansfile{trl,1};
                answer = cell2mat(answer);
                if answer == 'L'
                    area([(trial_len*num_trial)+(trial_len*(trl-1)) (trial_len*num_trial)+(trial_len*trl)], [2.5 2.5], 'FaceColor', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5)
                    hold on
                elseif answer == 'R'
                    area([(trial_len*num_trial)+(trial_len*(trl-1)) (trial_len*num_trial)+(trial_len*trl)], [1.5 1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5)
                end
            end
            area([0 trial_len*num_trial], [0.5 0.5], 'FaceColor', [0.6 0.2 0.8],  'EdgeColor', 'none', 'FaceAlpha', .5) % BR
             
            try % phys
                for i = 1:height(vert_triad_phys)
                    if rem(i,2) == 1; area([0.001*vert_triad_phys(i,1) 0.001*vert_triad_phys(i+1,1)], [-2.5 -2.5], 'FaceColor', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5); timestorage_phys_vert = timestorage_phys_vert + (vert_triad_phys(i+1,1) - vert_triad_phys(i,1)); end
                end
            catch
                area([0.001*vert_triad_phys(i,1) trial_len*num_trial], [-2.5 -2.5], 'FaceColor', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5);
                timestorage_phys_vert = timestorage_phys_vert + (1000*trial_len*num_trial - vert_triad_phys(i,1));
            end
            try
                for i = 1:height(horz_triad_phys)
                    if rem(i,2) == 1; area([0.001*horz_triad_phys(i,1) 0.001*horz_triad_phys(i+1,1)], [-1.5 -1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5); timestorage_phys_horz = timestorage_phys_horz + (horz_triad_phys(i+1,1) - horz_triad_phys(i,1)); end
                end
            catch
                area([0.001*horz_triad_phys(i,1) trial_len*num_trial], [-1.5 -1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5);
                timestorage_phys_horz = timestorage_phys_horz + (1000*trial_len*num_trial - horz_triad_phys(i,1));
            end
            try % BR
            for i = 1:height(vert_triad_br)
                if rem(i,2) == 1; area([(trial_len*num_trial)+0.001*vert_triad_br(i,1) (trial_len*num_trial)+0.001*vert_triad_br(i+1,1)], [-2.5 -2.5], 'FaceColor', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5); timestorage_br_vert = timestorage_br_vert + (vert_triad_br(i+1,1) - vert_triad_br(i,1)); end
            end
            catch
                area([(trial_len*num_trial)+0.001*vert_triad_br(i,1) (trial_len*num_trial)+trial_len*num_trial], [-2.5 -2.5], 'color', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5);
                timestorage_br_vert = timestorage_br_vert + (1000*trial_len*num_trial - vert_triad_br(i,1));
            end
            try
                for i = 1:height(horz_triad_br)
                    if rem(i,2) == 1; area([(trial_len*num_trial)+0.001*horz_triad_br(i,1) (trial_len*num_trial)+0.001*horz_triad_br(i+1,1)], [-1.5 -1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5); timestorage_br_horz = timestorage_br_horz + (horz_triad_br(i+1,1) - horz_triad_br(i,1)); end
                end
            catch
                area([(trial_len*num_trial)+0.001*horz_triad_br(i,1) (trial_len*num_trial)+trial_len*num_trial], [-1.5 -1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5);
                timestorage_br_horz = timestorage_br_horz + (1000*trial_len*num_trial - horz_triad_br(i,1));
            end
            rate_br_vert = timestorage_br_vert / (1000*trial_len*num_trial);
            rate_br_horz = timestorage_br_horz / (1000*trial_len*num_trial);
            rate_br_mix = 1 - rate_br_vert - rate_br_horz; if rate_br_mix < 0; rate_br_mix = 0; end
            rate_phys_vert = timestorage_phys_vert / (1000*trial_len*num_trial);
            rate_phys_horz = timestorage_phys_horz / (1000*trial_len*num_trial);
            ratestorage_br_vert = ratestorage_br_vert + rate_br_vert;
            ratestorage_br_horz = ratestorage_br_horz + rate_br_horz;
            ratestorage_br_mix = ratestorage_br_mix + rate_br_mix;
            ratestorage_phys_vert = ratestorage_phys_vert + rate_phys_vert;
            ratestorage_phys_horz = ratestorage_phys_horz + rate_phys_horz;
            
            xlabel('Time [s]')
            xticks(0:trial_len:trial_len*num_trial*2)
            xlim([0 trial_len*num_trial*2])
            ylabel('Report by button press         Stimulus switch')
            ylim([-3 3])
            yticks([-2.5 -1.5 0 0.5 1.5 2.5])
            yticklabels({'Red', 'Blue', ' ', 'BR','Blue','Red'})
            title({['superblock: ' num2str(spb) ', triads: ' num2str(trd)]; ...
                ['[BR] Red: ' num2str(100*rate_br_vert,3) ' %']; ...
                ['[BR] Blue: ' num2str(100*rate_br_horz,3) ' %']; ...
                ['[BR] Mixture: ' num2str(100*rate_br_mix,3) ' %']; ...
                ['[Phys] Red: ' num2str(100*rate_phys_vert,3) ' %']; ...
                ['[Phys] Blue: ' num2str(100*rate_phys_horz,3) ' %']});
            grid on
            
            % plot eye-traking
            if eye_track == 1
                yyaxis right
                % phys
                time = [];
                eyepos_x = [];
                eyepos_y = [];
                fix_loc = [];
                distance = [];
                eyeposfile_phys = readtable([phys_dir '/eyepos_' num2str(spb) '_' num2str(trd) '.csv']);
                fixlocfile_phys = readtable([phys_dir '/fixloc_' num2str(spb) '_' num2str(trd) '.csv']);
                for i = 1:height(eyeposfile_phys); time = vertcat(time, eyeposfile_phys{i,1}/1000); end
                for i = 1:height(eyeposfile_phys); eyepos_x = vertcat(eyepos_x, eyeposfile_phys{i,3}); end %deg
                for i = 1:height(eyeposfile_phys); eyepos_y = vertcat(eyepos_y, eyeposfile_phys{i,4}); end
                for i = 1:height(fixlocfile_phys); fix_loc(i,1) = fixlocfile_phys{i,1}; end
                for i = 1:height(fixlocfile_phys); fix_loc(i,2) = fixlocfile_phys{i,2}; end
                for i = 1:height(eyeposfile_phys)
                    fix_loc_l(i,1) = fix_loc(i,1);
                    if fix_loc_l(i,1) == 1 % fix_left/upleft
                        fix_loc_lx = fp_loc(1,1)+(fp_loc(1,3) - fp_loc(1,1))/2; % px
                        fix_loc_ly = fp_loc(1,2)+(fp_loc(1,4) - fp_loc(1,2))/2;
                    elseif fix_loc_l(i,1) == 2 % fix_right/upright
                        fix_loc_lx = fp_loc(2,1)+(fp_loc(2,3) - fp_loc(2,1))/2;
                        fix_loc_ly = fp_loc(2,2)+(fp_loc(2,4) - fp_loc(2,2))/2;
                    elseif fix_loc_l(i,1) == 3 % fix_up/belowleft
                        fix_loc_lx = fp_loc(3,1)+(fp_loc(3,3) - fp_loc(3,1))/2;
                        fix_loc_ly = fp_loc(3,2)+(fp_loc(3,4) - fp_loc(3,2))/2;
                    elseif fix_loc_l(i,1) == 4 % fix_below/belowright
                        fix_loc_lx = fp_loc(4,1)+(fp_loc(4,3) - fp_loc(4,1))/2;
                        fix_loc_ly = fp_loc(4,2)+(fp_loc(4,4) - fp_loc(4,2))/2;
                    end
                    fix_loc_r(i,1) = fix_loc(i,2);
                    if only_grating == 0
                        if fix_loc_r(i,1) == 1 % fix_left/upleft
                            fix_loc_rx = fp_loc(1,1)+(fp_loc(1,3) - fp_loc(1,1))/2; % px
                            fix_loc_ry = fp_loc(1,2)+(fp_loc(1,4) - fp_loc(1,2))/2;
                        elseif fix_loc_r(i,1) == 2 % fix_right/upright
                            fix_loc_rx = fp_loc(2,1)+(fp_loc(2,3) - fp_loc(2,1))/2;
                            fix_loc_ry = fp_loc(2,2)+(fp_loc(2,4) - fp_loc(2,2))/2;
                        elseif fix_loc_r(i,1) == 3 % fix_up/belowleft
                            fix_loc_rx = fp_loc(3,1)+(fp_loc(3,3) - fp_loc(3,1))/2;
                            fix_loc_ry = fp_loc(3,2)+(fp_loc(3,4) - fp_loc(3,2))/2;
                        elseif fix_loc_r(i,1) == 4 % fix_below/belowright
                            fix_loc_rx = fp_loc(4,1)+(fp_loc(4,3) - fp_loc(4,1))/2;
                            fix_loc_ry = fp_loc(4,2)+(fp_loc(4,4) - fp_loc(4,2))/2;
                        end
                    elseif only_grating == 1
                        if fix_loc_r(i,1) == 1 % fix_left/upleft
                            fix_loc_rx = fp_loc_grat(1,1)+(fp_loc_grat(1,3) - fp_loc_grat(1,1))/2; % px
                            fix_loc_ry = fp_loc_grat(1,2)+(fp_loc_grat(1,4) - fp_loc_grat(1,2))/2;
                        elseif fix_loc_r(i,1) == 2 % fix_right/upright
                            fix_loc_rx = fp_loc_grat(2,1)+(fp_loc_grat(2,3) - fp_loc_grat(2,1))/2;
                            fix_loc_ry = fp_loc_grat(2,2)+(fp_loc_grat(2,4) - fp_loc_grat(2,2))/2;
                        elseif fix_loc_r(i,1) == 3 % fix_up/belowleft
                            fix_loc_rx = fp_loc_grat(3,1)+(fp_loc_grat(3,3) - fp_loc_grat(3,1))/2;
                            fix_loc_ry = fp_loc_grat(3,2)+(fp_loc_grat(3,4) - fp_loc_grat(3,2))/2;
                        elseif fix_loc_r(i,1) == 4 % fix_below/belowright
                            fix_loc_rx = fp_loc_grat(4,1)+(fp_loc_grat(4,3) - fp_loc_grat(4,1))/2;
                            fix_loc_ry = fp_loc_grat(4,2)+(fp_loc_grat(4,4) - fp_loc_grat(4,2))/2;
                        end
                    end
                    % from px to deg|fixation spot
                    fix_loc_lx = fix_loc_lx - monitor_x/2; % px from centre
                    fix_loc_ly = fix_loc_ly - monitor_y/2;
                    fix_loc_rx = fix_loc_rx - monitor_x/2;
                    fix_loc_ry = fix_loc_ry - monitor_y/2;
                    width_cmperpx = screen_width/monitor_x; % cm/px of the screen
                    fix_loc_lx = fix_loc_lx * width_cmperpx; % cx from centre
                    fix_loc_ly = fix_loc_ly * width_cmperpx;
                    fix_loc_rx = fix_loc_rx * width_cmperpx;
                    fix_loc_ry = fix_loc_ry * width_cmperpx;
                    % V(deg) = 2*arctan(size/(2*DistFromMonitor))
                    fix_loc_lx = 2*atand(fix_loc_lx/(2*dist_scr)); % deg
                    fix_loc_ly = 2*atand(fix_loc_ly/(2*dist_scr));
                    fix_loc_rx = 2*atand(fix_loc_rx/(2*dist_scr));
                    fix_loc_ry = 2*atand(fix_loc_ry/(2*dist_scr));
                    
                    distance(i,1) = sqrt((eyepos_x(i,1)-fix_loc_lx).^2+(eyepos_y(i,1)-fix_loc_ly).^2); % L
                    distance(i,2) = sqrt((eyepos_x(i,1)-fix_loc_rx).^2+(eyepos_y(i,1)-fix_loc_ry).^2); % R
                end
                

                plot(time, distance(:,1), 'Color', 'red', 'LineStyle', '-');
                hold on
                plot(time, distance(:,2), 'Color', 'blue', 'LineStyle', '-');
                hold on
                
                time = [];
                eyepos_x = [];
                eyepos_y = [];
                fix_loc = [];
                distance = [];
                eyeposfile_bino = readtable([bino_dir '/eyepos_' num2str(spb) '_' num2str(trd) '.csv']);
                fixlocfile_bino = readtable([bino_dir '/fixloc_' num2str(spb) '_' num2str(trd) '.csv']);
                % bino
                for i = 1:height(eyeposfile_bino); time = vertcat(time, num_trial*trial_len+eyeposfile_bino{i,1}/1000); end
                for i = 1:height(eyeposfile_bino); eyepos_x = vertcat(eyepos_x, eyeposfile_bino{i,3}); end %deg
                for i = 1:height(eyeposfile_bino); eyepos_y = vertcat(eyepos_y, eyeposfile_bino{i,4}); end
                for i = 1:height(fixlocfile_bino); fix_loc(i,1) = fixlocfile_bino{i,1}; end
                for i = 1:height(fixlocfile_bino); fix_loc(i,2) = fixlocfile_bino{i,2}; end
                for i = 1:height(eyeposfile_bino)
                    fix_loc_l(i,1) = fix_loc(i,1);
                    if fix_loc_l(i,1) == 1 % fix_left/upleft
                        fix_loc_lx = fp_loc(1,1)+(fp_loc(1,3) - fp_loc(1,1))/2; % px
                        fix_loc_ly = fp_loc(1,2)+(fp_loc(1,4) - fp_loc(1,2))/2;
                    elseif fix_loc_l(i,1) == 2 % fix_right/upright
                        fix_loc_lx = fp_loc(2,1)+(fp_loc(2,3) - fp_loc(2,1))/2;
                        fix_loc_ly = fp_loc(2,2)+(fp_loc(2,4) - fp_loc(2,2))/2;
                    elseif fix_loc_l(i,1) == 3 % fix_up/belowleft
                        fix_loc_lx = fp_loc(3,1)+(fp_loc(3,3) - fp_loc(3,1))/2;
                        fix_loc_ly = fp_loc(3,2)+(fp_loc(3,4) - fp_loc(3,2))/2;
                    elseif fix_loc_l(i,1) == 4 % fix_below/belowright
                        fix_loc_lx = fp_loc(4,1)+(fp_loc(4,3) - fp_loc(4,1))/2;
                        fix_loc_ly = fp_loc(4,2)+(fp_loc(4,4) - fp_loc(4,2))/2;
                    end
                    fix_loc_r(i,1) = fix_loc(i,2);
                    if only_grating == 0
                        if fix_loc_r(i,1) == 1 % fix_left/upleft
                            fix_loc_rx = fp_loc(1,1)+(fp_loc(1,3) - fp_loc(1,1))/2; % px
                            fix_loc_ry = fp_loc(1,2)+(fp_loc(1,4) - fp_loc(1,2))/2;
                        elseif fix_loc_r(i,1) == 2 % fix_right/upright
                            fix_loc_rx = fp_loc(2,1)+(fp_loc(2,3) - fp_loc(2,1))/2;
                            fix_loc_ry = fp_loc(2,2)+(fp_loc(2,4) - fp_loc(2,2))/2;
                        elseif fix_loc_r(i,1) == 3 % fix_up/belowleft
                            fix_loc_rx = fp_loc(3,1)+(fp_loc(3,3) - fp_loc(3,1))/2;
                            fix_loc_ry = fp_loc(3,2)+(fp_loc(3,4) - fp_loc(3,2))/2;
                        elseif fix_loc_r(i,1) == 4 % fix_below/belowright
                            fix_loc_rx = fp_loc(4,1)+(fp_loc(4,3) - fp_loc(4,1))/2;
                            fix_loc_ry = fp_loc(4,2)+(fp_loc(4,4) - fp_loc(4,2))/2;
                        end
                    elseif only_grating == 1
                        if fix_loc_r(i,1) == 1 % fix_left/upleft
                            fix_loc_rx = fp_loc_grat(1,1)+(fp_loc_grat(1,3) - fp_loc_grat(1,1))/2; % px
                            fix_loc_ry = fp_loc_grat(1,2)+(fp_loc_grat(1,4) - fp_loc_grat(1,2))/2;
                        elseif fix_loc_r(i,1) == 2 % fix_right/upright
                            fix_loc_rx = fp_loc_grat(2,1)+(fp_loc_grat(2,3) - fp_loc_grat(2,1))/2;
                            fix_loc_ry = fp_loc_grat(2,2)+(fp_loc_grat(2,4) - fp_loc_grat(2,2))/2;
                        elseif fix_loc_r(i,1) == 3 % fix_up/belowleft
                            fix_loc_rx = fp_loc_grat(3,1)+(fp_loc_grat(3,3) - fp_loc_grat(3,1))/2;
                            fix_loc_ry = fp_loc_grat(3,2)+(fp_loc_grat(3,4) - fp_loc_grat(3,2))/2;
                        elseif fix_loc_r(i,1) == 4 % fix_below/belowright
                            fix_loc_rx = fp_loc_grat(4,1)+(fp_loc_grat(4,3) - fp_loc_grat(4,1))/2;
                            fix_loc_ry = fp_loc_grat(4,2)+(fp_loc_grat(4,4) - fp_loc_grat(4,2))/2;
                        end
                    end
                    % from px to deg|fixation spot
                    fix_loc_lx = fix_loc_lx - monitor_x/2; % px from centre
                    fix_loc_ly = fix_loc_ly - monitor_y/2;
                    fix_loc_rx = fix_loc_rx - monitor_x/2;
                    fix_loc_ry = fix_loc_ry - monitor_y/2;
                    width_cmperpx = screen_width/monitor_x; % cm/px of the screen
                    fix_loc_lx = fix_loc_lx * width_cmperpx; % cx from centre
                    fix_loc_ly = fix_loc_ly * width_cmperpx;
                    fix_loc_rx = fix_loc_rx * width_cmperpx;
                    fix_loc_ry = fix_loc_ry * width_cmperpx;
                    % V(deg) = 2*arctan(size/(2*DistFromMonitor))
                    fix_loc_lx = 2*atand(fix_loc_lx/(2*dist_scr)); % deg
                    fix_loc_ly = 2*atand(fix_loc_ly/(2*dist_scr));
                    fix_loc_rx = 2*atand(fix_loc_rx/(2*dist_scr));
                    fix_loc_ry = 2*atand(fix_loc_ry/(2*dist_scr));
                    
                    distance(i,1) = sqrt((eyepos_x(i,1)-fix_loc_lx).^2+(eyepos_y(i,1)-fix_loc_ly).^2); % L
                    distance(i,2) = sqrt((eyepos_x(i,1)-fix_loc_rx).^2+(eyepos_y(i,1)-fix_loc_ry).^2); % R
                end
                
                plot(time, distance(:,1), 'Color', 'red', 'LineStyle', '-');
                hold on
                plot(time, distance(:,2), 'Color', 'blue', 'LineStyle', '-');
                
                ylabel('Distance from fixation spot (deg)')
                ax1.YColor = 'k';
                ylim([0 2.5])
                legend('Dist from right eye spot', 'Dist from left eye spot')
            end
            
            filename = [fig_dir '/reco_' num2str(spb) '_' num2str(trd) '.png'];
            saveas(gcf,filename)
            filename = [fig_dir '/reco_' num2str(spb) '_' num2str(trd) '.fig'];
            saveas(gcf,filename)
            clf

       %% bino->phys
        elseif switchfile{trd,1} == 1 % if 0 phys->bino, if 1 bino->phys  
            t = tiledlayout(1,1);
            ax1 = axes(t);
            if eye_track == 1; yyaxis left; end
            % plot button press
            area([0 trial_len*num_trial], [0.5 0.5], 'FaceColor', [0.6 0.2 0.8],  'EdgeColor', 'none', 'FaceAlpha', .5) % BR
            hold on
            for trl = 1:num_trial
                answer = ansfile{trl,1};
                answer = cell2mat(answer);
                if answer == 'L'
                    area([(trial_len*num_trial)+(trial_len*(trl-1)) (trial_len*num_trial)+(trial_len*trl)], [2.5 2.5], 'FaceColor', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5)
                elseif answer == 'R'
                    area([(trial_len*num_trial)+(trial_len*(trl-1)) (trial_len*num_trial)+(trial_len*trl)], [1.5 1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5)
                end
            end
            
            
            try % BR
                for i = 1:height(vert_triad_br)
                    if rem(i,2) == 1; area([0.001*vert_triad_br(i,1) 0.001*vert_triad_br(i+1,1)], [-2.5 -2.5], 'FaceColor', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5); timestorage_br_vert = timestorage_br_vert + (vert_triad_br(i+1,1) - vert_triad_br(i,1)); end
                end
            catch
                area([0.001*vert_triad_br(i,1) trial_len*num_trial], [-2.5 -2.5], 'FaceColor', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5);
                timestorage_br_vert = timestorage_br_vert + (1000*trial_len*num_trial - vert_triad_br(i,1));
            end
            try
                for i = 1:height(horz_triad_br)
                    if rem(i,2) == 1; area([0.001*horz_triad_br(i,1) 0.001*horz_triad_br(i+1,1)], [-1.5 -1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5); timestorage_br_horz = timestorage_br_horz + (horz_triad_br(i+1,1) - horz_triad_br(i,1)); end
                end
            catch
                area([0.001*horz_triad_br(i,1) trial_len*num_trial], [-1.5 -1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5);
                timestorage_br_horz = timestorage_br_horz + (1000*trial_len*num_trial - horz_triad_br(i,1));
            end
            try % phys
                for i = 1:height(vert_triad_phys)
                    if rem(i,2) == 1; area([(trial_len*num_trial)+0.001*vert_triad_phys(i,1) (trial_len*num_trial)+0.001*vert_triad_phys(i+1,1)], [-2.5 -2.5], 'FaceColor', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5); timestorage_phys_vert = timestorage_phys_vert + (vert_triad_phys(i+1,1) - vert_triad_phys(i,1)); end
                end
            catch
                area([(trial_len*num_trial)+0.001*vert_triad_phys(i,1) (trial_len*num_trial)+trial_len*num_trial], [-2.5 -2.5], 'FaceColor', [0.8 0.2 0], 'EdgeColor', 'none', 'FaceAlpha', .5);
                timestorage_phys_vert = timestorage_phys_vert + (1000*trial_len*num_trial - vert_triad_phys(i,1));
            end
            try
                for i = 1:height(horz_triad_phys)
                    if rem(i,2) == 1; area([(trial_len*num_trial)+0.001*horz_triad_phys(i,1) (trial_len*num_trial)+0.001*horz_triad_phys(i+1,1)], [-1.5 -1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5); timestorage_phys_horz = timestorage_phys_horz + (horz_triad_phys(i+1,1) - horz_triad_phys(i,1)); end
                end
            catch
                area([(trial_len*num_trial)+0.001*horz_triad_phys(i,1) (trial_len*num_trial)+trial_len*num_trial], [-1.5 -1.5], 'FaceColor', [0 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', .5);
                timestorage_phys_horz = timestorage_phys_horz + (1000*trial_len*num_trial - horz_triad_phys(i,1));
            end
            rate_br_vert = timestorage_br_vert / (1000*trial_len*num_trial);
            rate_br_horz = timestorage_br_horz / (1000*trial_len*num_trial);
            rate_br_mix = 1 - rate_br_vert - rate_br_horz; if rate_br_mix < 0; rate_br_mix = 0; end
            rate_phys_vert = timestorage_phys_vert / (1000*trial_len*num_trial);
            rate_phys_horz = timestorage_phys_horz / (1000*trial_len*num_trial);
            ratestorage_br_vert = ratestorage_br_vert + rate_br_vert;
            ratestorage_br_horz = ratestorage_br_horz + rate_br_horz;
            ratestorage_br_mix = ratestorage_br_mix + rate_br_mix;
            ratestorage_phys_vert = ratestorage_phys_vert + rate_phys_vert;
            ratestorage_phys_horz = ratestorage_phys_horz + rate_phys_horz;
            
            xlabel('Time [s]')
            xticks(0:trial_len:trial_len*num_trial*2)
            xlim([0 trial_len*num_trial*2])
            ylabel('Report by button press         Stimulus switch')
            ax1.YColor = 'k';
            ylim([-3 3])
            yticks([-2.5 -1.5 0 0.5 1.5 2.5])
            yticklabels({'Red', 'Blue', ' ', 'BR','Blue','Red'})
            title({['superblock: ' num2str(spb) ', triads: ' num2str(trd)]; ...
                ['[BR] Red: ' num2str(100*rate_br_vert,3) ' %']; ...
                ['[BR] Blue: ' num2str(100*rate_br_horz,3) ' %']; ...
                ['[BR] Mixture: ' num2str(100*rate_br_mix,3) ' %']; ...
                ['[Phys] Red: ' num2str(100*rate_phys_vert,3) ' %']; ...
                ['[Phys] Blue: ' num2str(100*rate_phys_horz,3) ' %']});
            grid on
            
            % plot eye-traking
            if eye_track == 1
                yyaxis right
                time = [];
                eyepos_x = [];
                eyepos_y = [];
                fix_loc = [];
                distance = [];
                eyeposfile_bino = readtable([bino_dir '/eyepos_' num2str(spb) '_' num2str(trd) '.csv']);
                fixlocfile_bino = readtable([bino_dir '/fixloc_' num2str(spb) '_' num2str(trd) '.csv']);
                % bino
                for i = 1:height(eyeposfile_bino); time = vertcat(time, eyeposfile_bino{i,1}/1000); end
                for i = 1:height(eyeposfile_bino); eyepos_x = vertcat(eyepos_x, eyeposfile_bino{i,3}); end %deg
                for i = 1:height(eyeposfile_bino); eyepos_y = vertcat(eyepos_y, eyeposfile_bino{i,4}); end
                for i = 1:height(fixlocfile_bino); fix_loc(i,1) = fixlocfile_bino{i,1}; end
                for i = 1:height(fixlocfile_bino); fix_loc(i,2) = fixlocfile_bino{i,2}; end
                for i = 1:height(eyeposfile_bino)
                    fix_loc_l(i,1) = fix_loc(i,1);
                    if fix_loc_l(i,1) == 1 % fix_left/upleft
                        fix_loc_lx = fp_loc(1,1)+(fp_loc(1,3) - fp_loc(1,1))/2; % px
                        fix_loc_ly = fp_loc(1,2)+(fp_loc(1,4) - fp_loc(1,2))/2;
                    elseif fix_loc_l(i,1) == 2 % fix_right/upright
                        fix_loc_lx = fp_loc(2,1)+(fp_loc(2,3) - fp_loc(2,1))/2;
                        fix_loc_ly = fp_loc(2,2)+(fp_loc(2,4) - fp_loc(2,2))/2;
                    elseif fix_loc_l(i,1) == 3 % fix_up/belowleft
                        fix_loc_lx = fp_loc(3,1)+(fp_loc(3,3) - fp_loc(3,1))/2;
                        fix_loc_ly = fp_loc(3,2)+(fp_loc(3,4) - fp_loc(3,2))/2;
                    elseif fix_loc_l(i,1) == 4 % fix_below/belowright
                        fix_loc_lx = fp_loc(4,1)+(fp_loc(4,3) - fp_loc(4,1))/2;
                        fix_loc_ly = fp_loc(4,2)+(fp_loc(4,4) - fp_loc(4,2))/2;
                    end
                    fix_loc_r(i,1) = fix_loc(i,2);
                    if only_grating == 0
                        if fix_loc_r(i,1) == 1 % fix_left/upleft
                            fix_loc_rx = fp_loc(1,1)+(fp_loc(1,3) - fp_loc(1,1))/2; % px
                            fix_loc_ry = fp_loc(1,2)+(fp_loc(1,4) - fp_loc(1,2))/2;
                        elseif fix_loc_r(i,1) == 2 % fix_right/upright
                            fix_loc_rx = fp_loc(2,1)+(fp_loc(2,3) - fp_loc(2,1))/2;
                            fix_loc_ry = fp_loc(2,2)+(fp_loc(2,4) - fp_loc(2,2))/2;
                        elseif fix_loc_r(i,1) == 3 % fix_up/belowleft
                            fix_loc_rx = fp_loc(3,1)+(fp_loc(3,3) - fp_loc(3,1))/2;
                            fix_loc_ry = fp_loc(3,2)+(fp_loc(3,4) - fp_loc(3,2))/2;
                        elseif fix_loc_r(i,1) == 4 % fix_below/belowright
                            fix_loc_rx = fp_loc(4,1)+(fp_loc(4,3) - fp_loc(4,1))/2;
                            fix_loc_ry = fp_loc(4,2)+(fp_loc(4,4) - fp_loc(4,2))/2;
                        end
                    elseif only_grating == 1
                        if fix_loc_r(i,1) == 1 % fix_left/upleft
                            fix_loc_rx = fp_loc_grat(1,1)+(fp_loc_grat(1,3) - fp_loc_grat(1,1))/2; % px
                            fix_loc_ry = fp_loc_grat(1,2)+(fp_loc_grat(1,4) - fp_loc_grat(1,2))/2;
                        elseif fix_loc_r(i,1) == 2 % fix_right/upright
                            fix_loc_rx = fp_loc_grat(2,1)+(fp_loc_grat(2,3) - fp_loc_grat(2,1))/2;
                            fix_loc_ry = fp_loc_grat(2,2)+(fp_loc_grat(2,4) - fp_loc_grat(2,2))/2;
                        elseif fix_loc_r(i,1) == 3 % fix_up/belowleft
                            fix_loc_rx = fp_loc_grat(3,1)+(fp_loc_grat(3,3) - fp_loc_grat(3,1))/2;
                            fix_loc_ry = fp_loc_grat(3,2)+(fp_loc_grat(3,4) - fp_loc_grat(3,2))/2;
                        elseif fix_loc_r(i,1) == 4 % fix_below/belowright
                            fix_loc_rx = fp_loc_grat(4,1)+(fp_loc_grat(4,3) - fp_loc_grat(4,1))/2;
                            fix_loc_ry = fp_loc_grat(4,2)+(fp_loc_grat(4,4) - fp_loc_grat(4,2))/2;
                        end
                    end
                
                    % from px to deg|fixation spot
                    fix_loc_lx = fix_loc_lx - monitor_x/2; % px from centre
                    fix_loc_ly = fix_loc_ly - monitor_y/2;
                    fix_loc_rx = fix_loc_rx - monitor_x/2;
                    fix_loc_ry = fix_loc_ry - monitor_y/2;
                    width_cmperpx = screen_width/monitor_x; % cm/px of the screen
                    fix_loc_lx = fix_loc_lx * width_cmperpx; % cm from centre
                    fix_loc_ly = fix_loc_ly * width_cmperpx;
                    fix_loc_rx = fix_loc_rx * width_cmperpx;
                    fix_loc_ry = fix_loc_ry * width_cmperpx;
                    % V(deg) = 2*arctan(size/(2*DistFromMonitor))
                    fix_loc_lx = 2*atand(fix_loc_lx/(2*dist_scr)); % deg
                    fix_loc_ly = 2*atand(fix_loc_ly/(2*dist_scr));
                    fix_loc_rx = 2*atand(fix_loc_rx/(2*dist_scr));
                    fix_loc_ry = 2*atand(fix_loc_ry/(2*dist_scr));
                    
                    distance(i,1) = sqrt((eyepos_x(i,1)-fix_loc_lx).^2+(eyepos_y(i,1)-fix_loc_ly).^2); % L
                    distance(i,2) = sqrt((eyepos_x(i,1)-fix_loc_rx).^2+(eyepos_y(i,1)-fix_loc_ry).^2); % R
                end
                
                eye_track_left_plot = plot(time, distance(:,1), 'Color', 'red', 'LineStyle', '-', 'DisplayName','Dist from right eye spot');
                hold on
                eye_track_right_plot = plot(time, distance(:,2), 'Color', 'blue', 'LineStyle', '-', 'DisplayName','Dist from left eye spot');
                hold on
                legend([eye_track_left_plot,eye_track_right_plot])
                
                
                % phys
                time = [];
                eyepos_x = [];
                eyepos_y = [];
                fix_loc = [];
                distance = [];
                eyeposfile_phys = readtable([phys_dir '/eyepos_' num2str(spb) '_' num2str(trd) '.csv']);
                fixlocfile_phys = readtable([phys_dir '/fixloc_' num2str(spb) '_' num2str(trd) '.csv']);
                for i = 1:height(eyeposfile_phys); time = vertcat(time, num_trial*trial_len+eyeposfile_phys{i,1}/1000); end
                for i = 1:height(eyeposfile_phys); eyepos_x = vertcat(eyepos_x, eyeposfile_phys{i,3}); end %deg
                for i = 1:height(eyeposfile_phys); eyepos_y = vertcat(eyepos_y, eyeposfile_phys{i,4}); end
                for i = 1:height(fixlocfile_phys); fix_loc(i,1) = fixlocfile_phys{i,1}; end
                for i = 1:height(eyeposfile_phys)
                    if only_grating == 0
                        if fix_loc(i,1) == 1 % fix_left/upleft
                            fix_loc_x = fp_loc(1,1)+(fp_loc(1,3) - fp_loc(1,1))/2; % px
                            fix_loc_y = fp_loc(1,2)+(fp_loc(1,4) - fp_loc(1,2))/2;
                        elseif fix_loc(i,1) == 2 % fix_right/upright
                            fix_loc_x = fp_loc(2,1)+(fp_loc(2,3) - fp_loc(2,1))/2;
                            fix_loc_y = fp_loc(2,2)+(fp_loc(2,4) - fp_loc(2,2))/2;
                        elseif fix_loc(i,1) == 3 % fix_up/belowleft
                            fix_loc_x = fp_loc(3,1)+(fp_loc(3,3) - fp_loc(3,1))/2;
                            fix_loc_y = fp_loc(3,2)+(fp_loc(3,4) - fp_loc(3,2))/2;
                        elseif fix_loc(i,1) == 4 % fix_below/belowright
                            fix_loc_x = fp_loc(4,1)+(fp_loc(4,3) - fp_loc(4,1))/2;
                            fix_loc_y = fp_loc(4,2)+(fp_loc(4,4) - fp_loc(4,2))/2;
                        end
                    elseif only_grating == 1
                        if fix_loc(i,1) == 1 % fix_left/upleft
                            fix_loc_x = fp_loc_grat(1,1)+(fp_loc_grat(1,3) - fp_loc_grat(1,1))/2; % px
                            fix_loc_y = fp_loc_grat(1,2)+(fp_loc_grat(1,4) - fp_loc_grat(1,2))/2;
                        elseif fix_loc(i,1) == 2 % fix_right/upright
                            fix_loc_x = fp_loc_grat(2,1)+(fp_loc_grat(2,3) - fp_loc_grat(2,1))/2;
                            fix_loc_y = fp_loc_grat(2,2)+(fp_loc_grat(2,4) - fp_loc_grat(2,2))/2;
                        elseif fix_loc(i,1) == 3 % fix_up/belowleft
                            fix_loc_x = fp_loc_grat(3,1)+(fp_loc_grat(3,3) - fp_loc_grat(3,1))/2;
                            fix_loc_y = fp_loc_grat(3,2)+(fp_loc_grat(3,4) - fp_loc_grat(3,2))/2;
                        elseif fix_loc(i,1) == 4 % fix_below/belowright
                            fix_loc_x = fp_loc_grat(4,1)+(fp_loc_grat(4,3) - fp_loc_grat(4,1))/2;
                            fix_loc_y = fp_loc_grat(4,2)+(fp_loc_grat(4,4) - fp_loc_grat(4,2))/2;
                        end
                    end
                    % from px to deg|fixation spot
                    fix_loc_x = fix_loc_x - monitor_x/2; % px from centre
                    fix_loc_y = fix_loc_y - monitor_y/2;
                    
                    fix_loc_x = fix_loc_x * width_cmperpx; % cm from centre
                    fix_loc_y = fix_loc_y * width_cmperpx;
                    % V(deg) = 2*arctan(size/(2*DistFromMonitor))
                    fix_loc_x = 2*atand(fix_loc_x/(2*dist_scr)); % deg
                    fix_loc_y = 2*atand(fix_loc_y/(2*dist_scr));
                    
                    distance(i,1) = sqrt((eyepos_x(i,1)-fix_loc_x).^2+(eyepos_y(i,1)-fix_loc_y).^2);
                end
                

                eye_track_phys_plot = plot(time, distance(:,1), 'Color', 'k', 'LineStyle', '-',  'DisplayName','Dist from eye spot during physical switches');
                
                
                ylabel('Distance from fixation spot (deg)')
                ax1.YColor = 'k';
                ylim([0 3])
                legend([eye_track_left_plot,eye_track_right_plot,eye_track_phys_plot], 'Location', 'southoutside')
            end
            
            filename = [fig_dir '/reco_' num2str(spb) '_' num2str(trd) '.png'];
            saveas(gcf,filename)
            filename = [fig_dir '/reco_' num2str(spb) '_' num2str(trd) '.fig'];
            saveas(gcf,filename)
            clf
        end
    end
end 

sprintf('[BR] Red: %1f %%', (100*ratestorage_br_vert/(num_superblock*num_triad)))
sprintf('[BR] Blue: %1f %%', (100*ratestorage_br_horz/(num_superblock*num_triad)))
sprintf('[BR] Mixture: %1f %%', (100*ratestorage_br_mix/(num_superblock*num_triad)))
sprintf('[Phys] Red: %1f %%', (100*ratestorage_phys_vert/(num_superblock*num_triad)))
sprintf('[Phys] Blue: %1f %%', (100*ratestorage_phys_horz/(num_superblock*num_triad)))

       