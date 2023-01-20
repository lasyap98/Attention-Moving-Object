if Params.sub.subject==1
    block = [1, 2, 0, 1, 2, 0]; %face, house
elseif Params.sub.subject==2
    block = [0, 1, 2, 0, 1, 2]; %both, face, house
elseif Params.sub.subject==3
    block = [2, 0, 1, 2, 0, 1];
end

global blocksPerCon trialNo

trialNo=Params.trialNum;
blocksPerCon = 2;

house.eyedata = Block.EyeData(block==2);
house.fliptimes = Block.flipTimes(block==2);

face.eyedata = Block.EyeData(block==1);
face.fliptimes = Block.flipTimes(block==1);

both.eyedata = Block.EyeData(block==0);
both.fliptimes = Block.flipTimes(block==0);

for b = 1:Params.info.To_Blo-Params.info.FroBlo+1
    % stimuli positions across all trials in each block
    for iTrial = 1:trialNo
        currT_stim{b}{iTrial} = pT{Block.t2Run{b}(iTrial)}; 
        currT_stim2{b}{iTrial} = pT2{Block.t2Run{b}(iTrial)};
    end  
end

%% organize eye data
[house.eye_x, house.eye_y, house.eye_a, house.fixT] = organize_eyedata(house);
[face.eye_x, face.eye_y, face.eye_a, face.fixT] = organize_eyedata(face);
[both.eye_x, both.eye_y, both.eye_a, both.fixT] = organize_eyedata(both);

%% organize stim locations by block type
[house.face_x, house.face_y, house.house_x, house.house_y] = organize_stimLoc(find(block==2),currT_stim, currT_stim2);
[face.face_x, face.face_y, face.house_x, face.house_y] = organize_stimLoc(find(block==1),currT_stim, currT_stim2);
[both.face_x, both.face_y, both.house_x, both.house_y] = organize_stimLoc(find(block==0),currT_stim, currT_stim2);

%% cleanup ET data = stim and eye arrays will be equal
[house] = cleanup_eyedata(house);
[face] = cleanup_eyedata(face);
[both] = cleanup_eyedata(both);

%% find max correlation for optimal lag
% [house] = maxCorr_maxLag(house);
% [face] = maxCorr_maxLag(face);
% [both] = maxCorr_maxLag(both);

run analyze_LagETdata

%% sliding window of correlation: eye v stim

[house] = slidingwin_corr(house, 1, lag);
[face] = slidingwin_corr(face, 1, lag);
[both] = slidingwin_corr(both, 1, lag);

%% temporal correlation: stim1(face) v stim 2(house)

[house] = s1s2_corr(house,1);
[face] = s1s2_corr(face, 1);
[both] = s1s2_corr(both, 1);

%% save the variables

%if we're not in the subject folder in the main data folder, change pwd before saving the analysis variables
if ~contains(pwd, Params.dir) 
    cd(Params.dir)
end

save(sprintf('newanalyses_%s.mat', Params.sub.info), 'lag', 'Titrate', 'face', 'house', 'both', 'Params')

 %% temporal correlation: stim1(face) v stim 2(house)

function blocktype = s1s2_corr(blocktype, win_size)
global trialNo blocksPerCon

refresh = 1/60;
win_frames = round(win_size/refresh);


for b = 1:blocksPerCon
    
    stim1_x = blocktype.face_x{b};
    stim2_x = blocktype.house_x{b};
    stim1_y = blocktype.face_y{b};
    stim2_y = blocktype.house_y{b};
    
    for t = 1:trialNo

        win_shifts_s1s2x = length(stim2_x{t})-win_frames+1; %stim1 and stim2 will have same length
        win_shifts_s1s2y = length(stim2_y{t})-win_frames+1; %stim1 and stim2 will have same length
        
        assert(length(stim2_x{t})==length(stim1_x{t}), 'Array sizes don''t match!')
        assert(length(stim2_y{t})==length(stim1_y{t}), 'Array sizes don''t match!')

        corr_s1s2x{t} = zeros([1, win_shifts_s1s2x]);
        corr_s1s2y{t} = zeros([1, win_shifts_s1s2y]);

        for w = 1:win_shifts_s1s2x        
            corr_s1s2x{t}(w) = corr2(stim1_x{t}(w:w+win_frames-1), stim2_x{t}(w:w+win_frames-1));
        end 
        for w = 1:win_shifts_s1s2y       
            corr_s1s2y{t}(w) = corr2(stim1_y{t}(w:w+win_frames-1), stim2_y{t}(w:w+win_frames-1));
        end  

    end
    
    blocktype.corr_win.both_x{b} = corr_s1s2x;
    blocktype.corr_win.both_y{b} = corr_s1s2y;
    
    blocktype.mean_bothcorr{b} = zeros([2, trialNo]);
    
    for t = 1:length(corr_s1s2x)
        blocktype.mean_bothcorr{b}(1,t) = mean(corr_s1s2x{t});
        blocktype.mean_bothcorr{b}(2,t) = mean(corr_s1s2y{t});
    end
    
end

blocktype.mean_fh_x = mean([blocktype.mean_bothcorr{1}(1,:) blocktype.mean_bothcorr{2}(1,:)]);
blocktype.mean_fh_y = mean([blocktype.mean_bothcorr{1}(2,:) blocktype.mean_bothcorr{2}(2,:)]);


disp(['avg face-house corr: ' num2str(mean([blocktype.mean_fh_x, blocktype.mean_fh_y]))])
disp('done!')
end

%% sliding window of correlation: eye v stim

function blocktype = slidingwin_corr(blocktype, win_size, lag)
% EACH trial and EACH stim has an optimal lag
global trialNo blocksPerCon

refresh = 1/60;
win_frames = round(win_size/refresh);

for b = 1:blocksPerCon
    
    eye_x = blocktype.eye_x{b};
    eye_y = blocktype.eye_y{b};
    
    stim1_x = blocktype.face_x{b};
    stim1_y = blocktype.face_y{b};
    
    stim2_x = blocktype.house_x{b};
    stim2_y = blocktype.house_y{b};
    
    for t = 1:trialNo

        %1-1 lineup of eye-x and stim-1 and stim-2
        x1_trial = eye_x{t}(lag.avgmaxlag.facehouse:end);
        x2_trial = eye_x{t}(lag.avgmaxlag.facehouse:end);
        s1x_trial = stim1_x{t}(1:end-lag.avgmaxlag.facehouse+1);
        s2x_trial = stim2_x{t}(1:end-lag.avgmaxlag.facehouse+1);

        y1_trial = eye_y{t}(lag.avgmaxlag.facehouse:end);
        y2_trial = eye_y{t}(lag.avgmaxlag.facehouse:end);
        s1y_trial = stim1_y{t}(1:end-lag.avgmaxlag.facehouse+1);
        s2y_trial = stim2_y{t}(1:end-lag.avgmaxlag.facehouse+1);

        win_shifts_s1x = length(s1x_trial)-win_frames+1;
        win_shifts_s2x = length(s2x_trial)-win_frames+1;
        win_shifts_s1y = length(s1y_trial)-win_frames+1;
        win_shifts_s2y = length(s2y_trial)-win_frames+1;

        corr_win.face.x{t} = zeros([1, win_shifts_s1x]);
        corr_win.face.y{t} = zeros([1, win_shifts_s1y]);
        corr_win.house.x{t} = zeros([1, win_shifts_s2x]);
        corr_win.house.y{t} = zeros([1, win_shifts_s2y]);

        for w = 1:win_shifts_s1x        
            corr_win.face.x{t}(w) = corr2(x1_trial(w:w+win_frames-1), s1x_trial(w:w+win_frames-1));
        end
        for w = 1:win_shifts_s2x        
            corr_win.house.x{t}(w) = corr2(x2_trial(w:w+win_frames-1), s2x_trial(w:w+win_frames-1));
        end
        for w = 1:win_shifts_s1y        
            corr_win.face.y{t}(w) = corr2(y1_trial(w:w+win_frames-1), s1y_trial(w:w+win_frames-1));
        end
        for w = 1:win_shifts_s2y        
            corr_win.house.y{t}(w) = corr2(y2_trial(w:w+win_frames-1), s2y_trial(w:w+win_frames-1));
        end

        corrAnalyses.x1{b}{t} = x1_trial; %eye-x values shifted by optimal lag for stim1
        corrAnalyses.x2{b}{t} = x2_trial; %eye-x values shifted by optimal lag for stim2

        corrAnalyses.y1{b}{t} = y1_trial; %eye-y values shifted by optimal lag for stim1
        corrAnalyses.y2{b}{t} = y2_trial; %eye-y values shifted by optimal lag for stim2

        corrAnalyses.s1x{b}{t} = s1x_trial;
        corrAnalyses.s1y{b}{t} = s1y_trial;

        corrAnalyses.s2x{b}{t} = s2x_trial;
        corrAnalyses.s2y{b}{t} = s2y_trial;

    end

    blocktype.corrwin{b} = corr_win;
    
    blocktype.mean_face{b} = zeros([2, trialNo]);
    blocktype.mean_house{b} = zeros([2, trialNo]);
    
    for t = 1:trialNo
        blocktype.mean_face{b}(1,t) = mean(corr_win.face.x{t});
        blocktype.mean_face{b}(2,t) = mean(corr_win.face.y{t});
        blocktype.mean_house{b}(1,t) = mean(corr_win.house.x{t});
        blocktype.mean_house{b}(2,t) = mean(corr_win.house.y{t});
    end
   
end

%the inputs for each sliding window correlation output, shifted by the lag
blocktype.corr_inputs = corrAnalyses;

%calculate the average correlation for each stimulus across blocks
blocktype.meanfacex = mean([blocktype.mean_face{1}(1,:) blocktype.mean_face{2}(1,:)]);
blocktype.meanfacey = mean([blocktype.mean_face{1}(2,:) blocktype.mean_face{2}(2,:)]);
blocktype.meanhousex = mean([blocktype.mean_house{1}(1,:) blocktype.mean_house{2}(1,:)]);
blocktype.meanhousey = mean([blocktype.mean_house{1}(2,:) blocktype.mean_house{2}(2,:)]);

%calculate the average correlation for each stimulus across x/y    
blocktype.meanfacexy = mean([blocktype.meanfacex, blocktype.meanfacey]);
blocktype.meanhousexy = mean([blocktype.meanhousex, blocktype.meanhousey]);

disp(['avg house corr: ' num2str(mean([blocktype.meanhousex, blocktype.meanhousey]))])
disp(['avg face corr: ' num2str(mean([blocktype.meanfacex, blocktype.meanfacey]))])
end

%% find max correlation for optimal lag - replace this with lag trial analysis
% eye- and stim- arrays are all the same size

function blocktype = maxCorr_maxLag(blocktype)

for b = 1:length(blocktype.eye_x)
    for t = 1:length(blocktype.eye_x{b})
        for n=1:30
            s1x_temp = blocktype.face_x{b}{t}(1:end+1-n);
            s1y_temp = blocktype.face_y{b}{t}(1:end+1-n);

            s2x_temp = blocktype.house_x{b}{t}(1:end+1-n);
            s2y_temp = blocktype.house_y{b}{t}(1:end+1-n);

            ex_temp = blocktype.eye_x{b}{t}(n:end);
            ey_temp = blocktype.eye_y{b}{t}(n:end);

            blocktype.r_score.face{b}{t}(n,1) = corr2(s1x_temp, ex_temp);
            blocktype.r_score.face{b}{t}(n,2) = corr2(s1y_temp, ey_temp);

            blocktype.r_score.house{b}{t}(n,1) = corr2(s2x_temp, ex_temp);
            blocktype.r_score.house{b}{t}(n,2) = corr2(s2y_temp, ey_temp);

        end

        blocktype.max_corr.face{b}(t,1) = max(blocktype.r_score.face{b}{t}(:,1)); % max correlation for x values
        blocktype.max_corr.face{b}(t,2) = max(blocktype.r_score.face{b}{t}(:,2)); % max correlation for y values

        blocktype.max_corr.house{b}(t,1) = max(blocktype.r_score.house{b}{t}(:,1)); % max correlation for x values
        blocktype.max_corr.house{b}(t,2) = max(blocktype.r_score.house{b}{t}(:,2)); % max correlation for y values

        blocktype.max_lag.face{b}(t,1) = find(blocktype.r_score.face{b}{t}(:,2)==max(blocktype.r_score.face{b}{t}(:,2))); % optimal lag for x 
        blocktype.max_lag.face{b}(t,2) = find(blocktype.r_score.face{b}{t}(:,1)==max(blocktype.r_score.face{b}{t}(:,1))); % optimal lag for y

        blocktype.max_lag.house{b}(t,1) = find(blocktype.r_score.house{b}{t}(:,2)==max(blocktype.r_score.house{b}{t}(:,2))); % optimal lag for x 
        blocktype.max_lag.house{b}(t,2) = find(blocktype.r_score.house{b}{t}(:,1)==max(blocktype.r_score.house{b}{t}(:,1))); % optimal lag for y

    end
end

end
%% remove frames that had no eye data

function blocktype = cleanup_eyedata(blocktype)

global trialNo blocksPerCon 

blink_buffer = 5;
blink_spikeSize = 50;

for b = 1:blocksPerCon
    for trial = 1:trialNo
        
        %closed pupil before fixation start?
        zero_pupil1 = find(blocktype.eye_a{b}{trial}==0);
        
        %how much extra ETdata was collected before second stim frame was presented (trial started)
        beforeFix = find(blocktype.fixT{b}{trial}<blocktype.fliptimes{b}{trial}(2)); 
        
        if blocktype.eye_a{b}{trial}(beforeFix(end-1))==0
            %for x-y we don't want to remove the same ones as eye-a in case there was 
            zero_pupilb4fix = intersect(zero_pupil1, beforeFix); %find when zero pupil before fixation 
            nozero_b4fix = beforeFix(1:(end-length(zero_pupilb4fix)+1)); %remove before-fixation points to account for no pupil frames, use for eye values
        else
            zero_pupilb4fix = intersect(zero_pupil1, beforeFix); %find when zero pupil before fixation 
            nozero_b4fix = beforeFix(1:(end-length(zero_pupilb4fix))); %remove before-fixation points to account for no pupil frames, use for eye values
        end
        
        blocktype.eye_x{b}{trial} = blocktype.eye_x{b}{trial}(nozero_b4fix(end-1):end);
        blocktype.eye_y{b}{trial} = blocktype.eye_y{b}{trial}(nozero_b4fix(end-1):end);
        blocktype.eye_a{b}{trial} = blocktype.eye_a{b}{trial}(beforeFix(end-1):end);
        
        %where do we not see the pupil? (after fixation now)
        zero_pupil2 = find(blocktype.eye_a{b}{trial}==0);
        
        
        del = 0; %keep track of how many we removing with a=0
        for frame = zero_pupil2
            blocktype.face_x{b}{trial}(frame-del) = [];
            blocktype.house_x{b}{trial}(frame-del) = [];
            blocktype.face_y{b}{trial}(frame-del) = [];
            blocktype.house_y{b}{trial}(frame-del) = [];
            del = del+1;
        end
        
        %remove disturbances from blinks (spikes in eye movement)
        %% with x values
        a = sort([find(diff(blocktype.eye_x{b}{trial})>blink_spikeSize) ...
            find(diff(blocktype.eye_x{b}{trial})<-blink_spikeSize)]);
        
        if a %if there are actually blinks
            del = 0;
            d_a = find(diff(a)>3);
            i_d = 1;
            i_a = 1;
            while i_d <= length(d_a)+1
                if i_d == length(d_a)+1
                    delete = a(i_a)-del-blink_buffer : a(end)-del+blink_buffer;
                    blocktype.eye_x{b}{trial}(delete>0) = [];
                    blocktype.face_x{b}{trial}(delete>0) = [];
                    blocktype.house_x{b}{trial}(delete>0) = [];
                    del = del+length(delete(delete>0));
                    i_d= i_d+1;
                else
                    delete = a(i_a)-del-blink_buffer: a(d_a(i_d))-del+blink_buffer;
                    blocktype.eye_x{b}{trial}(delete>0) = [];
                    blocktype.face_x{b}{trial}(delete>0) = [];
                    blocktype.house_x{b}{trial}(delete>0) = [];
                    del = del+length(delete(delete>0));
                    i_a = d_a(i_d)+1;
                    i_d= i_d+1;
                end
            end
        end
        
        %% with y values
        a = sort([find(diff(blocktype.eye_y{b}{trial})>blink_spikeSize)...
            find(diff(blocktype.eye_y{b}{trial})<-blink_spikeSize)]);
        %if there are actually blinks
        if a
            del = 0;
            d_a = find(diff(a)>3);
            i_d = 1;
            i_a = 1;
            while i_d <= length(d_a)+1
                if i_d == length(d_a)+1
                    delete = a(i_a)-del-blink_buffer : a(end)-del+blink_buffer;
                    blocktype.eye_y{b}{trial}(delete>0) = [];
                    blocktype.face_y{b}{trial}(delete>0) = [];
                    blocktype.house_y{b}{trial}(delete>0) = [];
                    del = del+length(delete(delete>0));
                    i_d= i_d+1;
                else
                    delete = a(i_a)-del-blink_buffer: a(d_a(i_d))-del+blink_buffer;
                    blocktype.eye_y{b}{trial}(delete>0) = [];
                    blocktype.face_y{b}{trial}(delete>0) = [];
                    blocktype.house_y{b}{trial}(delete>0) = [];
                    del = del+length(delete(delete>0));
                    i_a = d_a(i_d)+1;
                    i_d= i_d+1;
                end
            end
        end
        %%
        assert(length(blocktype.face_x{b}{trial})==length(blocktype.eye_x{b}{trial}), 'Array sizes don''t match!')
        assert(length(blocktype.face_y{b}{trial})==length(blocktype.eye_y{b}{trial}), 'Array sizes don''t match!')

    end
end

disp('all cleaned up!')
end
%% organized stim coordinates

function [stim1_x, stim1_y, stim2_x, stim2_y] = organize_stimLoc(whichBlocks, currT_stim, currT_stim2)
% stim1 = face
% stim2 = house
global blocksPerCon

for b = 1:blocksPerCon
    % stimulus1-x values across time
    stim1_x{b} = {};
    for trial = 1:length(currT_stim{whichBlocks(b)})
        stim1_x{b}{trial} = [];
        for i = 1:length(currT_stim{whichBlocks(b)}{trial})
            stim1_x{b}{trial}(end+1)= mean([currT_stim{whichBlocks(b)}{trial}(i,1),currT_stim{whichBlocks(b)}{trial}(i,3)]);
        end
    end
end

for b = 1:blocksPerCon
    % stimulus1-y values across time
    stim1_y{b} = {};
    for trial = 1:length(currT_stim{whichBlocks(b)})
        stim1_y{b}{trial} = [];
        for i = 1:length(currT_stim{whichBlocks(b)}{trial})
            stim1_y{b}{trial}(end+1)= mean([currT_stim{whichBlocks(b)}{trial}(i,2),currT_stim{whichBlocks(b)}{trial}(i,4)]);
        end
    end
end

for b = 1:blocksPerCon
    % stimulus2-x values across time
    stim2_x{b} = {};
    for trial = 1:length(currT_stim2{whichBlocks(b)})
        stim2_x{b}{trial} = [];
        for i = 1:length(currT_stim2{whichBlocks(b)}{trial})
            stim2_x{b}{trial}(end+1)= mean([currT_stim2{whichBlocks(b)}{trial}(i,1),currT_stim2{whichBlocks(b)}{trial}(i,3)]);
        end
    end
end

for b = 1:blocksPerCon
    % stimulus2-y values across time
    stim2_y{b} = {};
    for trial = 1:length(currT_stim2{whichBlocks(b)})
        stim2_y{b}{trial} = [];
        for i = 1:length(currT_stim2{whichBlocks(b)}{trial})
            stim2_y{b}{trial}(end+1)= mean([currT_stim2{whichBlocks(b)}{trial}(i,2),currT_stim2{whichBlocks(b)}{trial}(i,4)]);
        end
    end
end
end

%% organize all eye x and y values by trial
function [eye_x,eye_y,eye_a,fix_t] = organize_eyedata(blocktype)

global trialNo blocksPerCon

for b = 1:blocksPerCon
    for trial = 1:trialNo
    % eyetracking-x values across time
        eye_x{b}{trial} = blocktype.eyedata{b}.mx{trial}; %30+ get values after fixating at center for 500ms%
    end
end
for b = 1:blocksPerCon
    for trial = 1:trialNo
    % eyetracking-y values across time
        eye_y{b}{trial} = blocktype.eyedata{b}.my{trial}; %30+ get values after fixating at center for 500ms%
    end
end
for b = 1:blocksPerCon
    for trial = 1:trialNo
    % eyetracking-a values across time
        eye_a{b}{trial} = blocktype.eyedata{b}.ma{trial}; %30+ get values after fixating at center for 500ms%
    end
end

for b = 1:blocksPerCon
    for trial = 1:trialNo
    % eyetracking-a values across time
        fix_t{b}{trial} = blocktype.eyedata{b}.FixDoneT{trial}; %30+ get values after fixating at center for 500ms%
    end
end
end
