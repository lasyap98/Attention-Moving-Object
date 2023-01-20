function corr_idx = calculate_stimCorr(possTrials, possTrials2)

global Params Display
%% calculate average temporal correlation between stim1 and stim2
% select the lowest correlation trajectories across all possible trajectories

% stimulus1-x values across time
for trial = 1:length(possTrials)
    stim1_x{trial} = [];
    for i = 1:length(possTrials{trial})
        stim1_x{trial}(end+1)= mean([possTrials{trial}(i,1),possTrials{trial}(i,3)]);
    end
end

% stimulus1-y values across time
for trial = 1:length(possTrials)
    stim1_y{trial} = [];
    for i = 1:length(possTrials{trial})
        stim1_y{trial}(end+1)= mean([possTrials{trial}(i,2),possTrials{trial}(i,4)]);
    end
end

% stimulus2-x values across time
for trial = 1:length(possTrials2)
    stim2_x{trial} = [];
    for i = 1:length(possTrials2{trial})
        stim2_x{trial}(end+1)= mean([possTrials2{trial}(i,1),possTrials2{trial}(i,3)]);
    end
end

% stimulus2-y values across time
for trial = 1:length(possTrials2)
    stim2_y{trial} = [];
    for i = 1:length(possTrials2{trial})
        stim2_y{trial}(end+1)= mean([possTrials2{trial}(i,2),possTrials2{trial}(i,4)]);
    end
end

%%

win_size = 1; %sec
win_frames = round(win_size/Display.refresh);

for t = 1:length(stim1_x)

    win_shifts_s1s2 = length(stim2_x{t})-win_frames+1; %stim1 and stim2 will have same length, x or y
    
    corr_win.s1s2.x{t} = zeros([1, win_shifts_s1s2]);
    corr_win.s1s2.y{t} = zeros([1, win_shifts_s1s2]);
    corr_win.s1s2.xy{t} = zeros([1, win_shifts_s1s2]);
    
    for w = 1:win_shifts_s1s2        
        corr_win.s1s2.x{t}(w) = corr2(stim1_x{t}(w:w+win_frames-1), stim2_x{t}(w:w+win_frames-1));
    end 
    for w = 1:win_shifts_s1s2        
        corr_win.s1s2.y{t}(w) = corr2(stim1_y{t}(w:w+win_frames-1), stim2_y{t}(w:w+win_frames-1));
    end  
end

mean_s1s2corr = zeros([2, Params.trialNum]); %find the average correlation for each possible trial
abs_mean_s1s2corr2 = zeros([1, Params.trialNum]); %find the average correlation for each possible trial
mean_s1s2corr2 = zeros([1, Params.trialNum]);

for t = 1:length(corr_win.s1s2.x)
    mean_s1s2corr(1,t) = mean(corr_win.s1s2.x{t});
    mean_s1s2corr(2,t) = mean(corr_win.s1s2.y{t});
    abs_mean_s1s2corr2(t) = mean([abs(mean(corr_win.s1s2.x{t})), abs(mean(corr_win.s1s2.y{t}))]); %average of absolute x and y correlations for each trial
    mean_s1s2corr2(t) = mean([corr_win.s1s2.x{t}, corr_win.s1s2.y{t}]);
end


%% select the lowest correlating trajectories

[abs_corrmeans_sorted,sort_i] = sort(abs_mean_s1s2corr2); %sort based on absolute distance from 0
% corrmeans_sorted = mean_s1s2corr2(sort_i); %sort the original corr averages the same way
lowestCorr_vals  = abs_corrmeans_sorted(1:Params.trialNum*2); %select the first n=Params.trialNum lowest corr trials
%maybe reference more because we'll have to remove invalid (out-of-bounds) trajectories

corr_idx = [];
for low = 1:Params.trialNum*2
    %corresponding idx with lowest corr trials
    corr_idx(low) = find(abs_mean_s1s2corr2==lowestCorr_vals(low)); %use this to reference possTrials
end
    
end