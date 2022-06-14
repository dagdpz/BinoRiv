% Plot sample intervals of eye-tracking data 
% Author: Ryo Segawa (whizznihil.kid@gmail.com)
function sample_mass = sample_intervals()
close all

subj_type = 0;
subj = 'Rebecca_0609';
num_superblock = 1;
num_triad = 3;

if subj_type == 0
    subj_dir = fullfile('recording/human/', subj);
elseif subj_type == 1
    subj_dir = fullfile('recording/monkey/', subj);
end
phys_dir = [subj_dir '/report/phys'];
bino_dir = [subj_dir '/report/bino'];

% eyeposfile = readtable([phys_dir '/eyepos_' num2str(num_superblock) '_' num2str(num_triad) '.csv']);
eyeposfile = readtable([bino_dir '/eyepos_' num2str(num_superblock) '_' num2str(num_triad) '.csv']);
% eyeposfile = readtable([phys_dir '/eyepos_1_2.csv']);
% eyeposfile = readtable([bino_dir '/eyepos_1_2.csv']);
% eyeposfile = readtable([bino_dir '/test.csv']);


time = [];
% for i = 1:height(eyeposfile); time = vertcat(time, num_trial*trial_len+eyeposfile{i,1}/1000); end
for i = 1:height(eyeposfile); time = vertcat(time, eyeposfile{i,1}/1000); end

sample_interval = [];
try
    for i = 1:length(time)
        sample_interval = vertcat(sample_interval, 1000*(time(i+1,1)-time(i,1)));
        sample_interval(i,:) = round(sample_interval(i,:), 2);
    end
catch
end

sample_interval_sort = sort(sample_interval);
array_repeated = arrayfun(@(x)length(find(sample_interval_sort == x)), unique(sample_interval_sort), 'Uniform', false);
array_repeated = cell2mat(array_repeated);
for i=1:length(sample_interval_sort)
    try
        while sample_interval_sort(i,:) == sample_interval_sort(i+1,:)
            sample_interval_sort(i+1,:)=[];
        end
    catch
    end
end

    
bar(sample_interval_sort, array_repeated)
xlim([0.5 1.5])
xlabel('Sampling interval (ms)')
ylabel('Num of samples')

[M,I] = max(array_repeated);
sample_mass = sample_interval_sort(I,:);