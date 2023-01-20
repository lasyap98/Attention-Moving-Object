% load IN THIS ORDER: 
% for lag, EndPractice 
% for main, main data >> newanalyses files
% clearvars -except face house both Params Block Lag block el edf_main edf_lag messages EL_lag face_EL flower_EL both_EL

edf_lag = Edf2Mat(sprintf('%s_1.edf', Params.sub.info)); %load the lag/practice data file
edf_main = Edf2Mat(sprintf('%s_2.edf', Params.sub.info)); %load the main expt data file

save(sprintf('edf_analyses_%s.mat', Params.sub.info), 'messages', ...
    'trial_msgs', 'trial_msgs_tm', 'EL_lag', 'face_EL', 'flower_EL', 'both_EL')
%% get the lag/prac messages

msgs_lag = edf_lag.Events.Messages.info;
msg_times_lag = edf_lag.Events.Messages.time;

%what messages did we pass on to Eyelink in the first place? - we denoted the start and end of each block, but not for each trial!
%to get each specific trial times, go to Eyelink.t data
% sprintf('EVENT_LagStart_%s', s_name)
% sprintf('EVENT_LagEnd_%s', s_name)
% sprintf('EVENT_PracRun%d_Start', Prac_run)
% sprintf('EVENT_PracRun%d_End', Prac_run)
% sprintf('EVENT_PracticeEnd')


e = 1;
for m = 1:length(msgs_lag)
    if contains(msgs_lag{m}, 'EVENT_Lag') %gather all lag event messages
        evt_msgs_lag{e} = msgs_lag{m};
        evt_msgs_lag_tm{e} = msg_times_lag(m); %timestamps of each message
        e = e+1;
    end
end

e = 1;
for m = 1:length(msgs_lag)
    if contains(msgs_lag{m}, 'EVENT_Prac') %gather all practice event messages 
        evt_msgs_prac{e} = msgs_lag{m};
        evt_msgs_prac_tm{e} = msg_times_lag(m); %timestamps of each message
        e = e+1;
    end
end

messages.msgs_lag = msgs_lag;
messages.msgs_times_lag = msg_times_lag;
messages.evt_msgs_lag = evt_msgs_lag;
messages.evt_msgs_lag_tm = evt_msgs_lag_tm;
messages.evt_msgs_prac = evt_msgs_prac;
messages.evt_msgs_prac_tm = evt_msgs_prac_tm;

% return msgs_lag, msgs_times_lag, evt_msgs_lag, evt_msgs_lag_tm, evt_msgs_prac, evt_msgs_prac_tm

% clean up eye data from Eyelink
for b = 1:s
    for t = 1:Lag.tperS
        % find where fixation ends and trial begins: 60 Hz
        afterFix_idx = find(Lag.EyeData{b}.FixDoneT{t} >= Lag.flipTimes{b}{t}(2)); %which frames are after fixation
        endOfFix_id = afterFix_idx(1); %the first frame after fixation
        
        % find the start/end timestamp idx in Samples timeline based on  60 Hz
        EL_lag.startTrial_abs_id{b}{t} = find(edf_lag.Samples.time==Lag.EyeData{b}.time{t}(1)); %trial start before fixation
        EL_lag.startTrial_id{b}{t} = find(edf_lag.Samples.time==Lag.EyeData{b}.time{t}(endOfFix_id)); %first point after fixation
        EL_lag.endTrial_id{b}{t} = find(edf_lag.Samples.time==Lag.EyeData{b}.time{t}(end)); %end of trial
        
        % the range of entire trial we Sample
        EL_lag.rangeTrial_id{b}{t} = EL_lag.startTrial_id{b}{t}:EL_lag.endTrial_id{b}{t};
        EL_lag.rangeTrial_t{b}{t} = edf_lag.Samples.time(EL_lag.rangeTrial_id{b}{t});
        EL_lag.trialDur_nocuts{b}{t}= length(EL_lag.rangeTrial_id{b}{t}); % in milliseconds, also the number of eye position points we start off with (not counting binks)
        
        % get all Eyelink values that are not missing data
        gx = edf_lag.Samples.gx(edf_lag.Samples.gx~=el.MISSING_DATA)';
        gy = edf_lag.Samples.gy(edf_lag.Samples.gy~=el.MISSING_DATA)';
        pa = edf_lag.Samples.pa(edf_lag.Samples.gy~=el.MISSING_DATA)';
        
        % get the values: keep a copy of the uncut values
        EL_lag.x_full{b}{t} = gx(EL_lag.rangeTrial_id{b}{t});
        EL_lag.y_full{b}{t} = gy(EL_lag.rangeTrial_id{b}{t});
        EL_lag.a_full{b}{t} = pa(EL_lag.rangeTrial_id{b}{t});
        EL_lag.t_full{b}{t} = edf_lag.Samples.time(EL_lag.rangeTrial_id{b}{t})';
        
        % make a copy of values to edit so we have 2 versions: edited (below) and unedited (above)
        EL_lag.x{b}{t} = EL_lag.x_full{b}{t};
        EL_lag.y{b}{t} = EL_lag.y_full{b}{t};
        EL_lag.a{b}{t} = EL_lag.a_full{b}{t};
        EL_lag.t{b}{t} = EL_lag.t_full{b}{t};
        
        % find the blinks
        
        % use Eyelink's built-in software to detect blinks, will remove the same
        % regions from x and y arrays
        EL_lag.blink_id{b}{t} = find(edf_lag.Events.Eblink.start>EL_lag.rangeTrial_t{b}{t}(1) ...
            & edf_lag.Events.Eblink.start<EL_lag.rangeTrial_t{b}{t}(end)); %is there a blink that starts during trial
        %the 'corrupted' eye data before and after a blink to remove (ms)
        blink_buffer = 100;
        EL_lag.to_remove{b}{t} = []; %add these indices to the ones we need to remove from the stimulus
        
        if ~isempty(EL_lag.blink_id{b}{t})
            for blink = EL_lag.blink_id{b}{t}
                blink_start_id = find(EL_lag.rangeTrial_t{b}{t}==edf_lag.Events.Eblink.start(blink)); % index of blink start time
                blink_start_t = edf_lag.Events.Eblink.start(blink); % blink start time
                blink_end_id = find(EL_lag.rangeTrial_t{b}{t}==edf_lag.Events.Eblink.end(blink)); %index of blink end time
                blink_end_t = edf_lag.Events.Eblink.end(blink); % blink end time
                
                % to make sure the buffer doesn't go beyond the trial bounds
                if blink_start_id-blink_buffer > 0
                    st_rem = blink_start_id-blink_buffer;
                else
                    st_rem = 1;
                end
                if blink_end_id+blink_buffer <= length(EL_lag.t{b}{t})
                    end_rem = blink_end_id+blink_buffer;
                else
                    end_rem = length(EL_lag.t{b}{t}); % x and y data will be same length
                end
                EL_lag.to_remove{b}{t} = [EL_lag.to_remove{b}{t} st_rem:end_rem]; %the 'corrupted' eye data before and after a blink to remove (ms)
                
                
            end
        end
        
        % deletion: do it only once per trial or displacement will be wrong
        EL_lag.x{b}{t}(EL_lag.to_remove{b}{t}) = [];
        EL_lag.y{b}{t}(EL_lag.to_remove{b}{t}) = [];
        EL_lag.a{b}{t}(EL_lag.to_remove{b}{t}) = [];
        EL_lag.t{b}{t}(EL_lag.to_remove{b}{t}) = [];
        
        % to check that all x/y nans were included in removal
        x_isnan = isnan(EL_lag.x_full{b}{t});
        assert(sum(x_isnan==isnan(EL_lag.y_full{b}{t}))==length(x_isnan), 'NaNs don''t match!')
        assert(sum(ismember(find(x_isnan), EL_lag.to_remove{b}{t}))==sum(x_isnan), 'NaNs not deleted!')
        
    end
end

% get the stimulus data


% calculate average stimulus x and y values (center points) - (s1:face, s2:flower)
for t = 1:Lag.tperS

    EL_lag.s1_x_600{t} = []; EL_lag.s1_y_600{t} = [];
    EL_lag.s2_x_600{t} = []; EL_lag.s2_y_600{t} = [];
    % ONLY gives us 600 points
    
    for i = 1:length(Lag.currTrial_stim2{t})
        EL_lag.s1_x_600{t}(end+1)= mean([Lag.currTrial_stim{t}(i,1),Lag.currTrial_stim{t}(i,3)]);
        EL_lag.s1_y_600{t}(end+1)= mean([Lag.currTrial_stim{t}(i,2),Lag.currTrial_stim{t}(i,4)]);
        EL_lag.s2_x_600{t}(end+1)= mean([Lag.currTrial_stim2{t}(i,1),Lag.currTrial_stim2{t}(i,3)]);
        EL_lag.s2_y_600{t}(end+1)= mean([Lag.currTrial_stim2{t}(i,2),Lag.currTrial_stim2{t}(i,4)]);
    end
end

% need to match the number of stimulus points to eye data points: linear interpolation
for b = 1:s
    for t = 1:Lag.tperS
        % x and y should have the same number of points
        xi = linspace(1,600, EL_lag.trialDur_nocuts{b}{t}); % 1 to 'length(stim)', with 'length of uncut trial' amount of points
        
        if b==1
%             EL_lag.s1_x{t} = []; EL_lag.s1_y{t} = [];
            
            %interpolate the data to full length of uncut trial
            EL_lag.s1_x{t} = interp1(1:600, EL_lag.s1_x_600{t}, xi);
            EL_lag.s1_y{t} = interp1(1:600, EL_lag.s1_y_600{t}, xi);
            
            %remove the eye data (nan + blinks)
            EL_lag.s1_x{t}(EL_lag.to_remove{b}{t}) = []; EL_lag.s1_y{t}(EL_lag.to_remove{b}{t}) = [];
            
            %make sure we have the same length arrays for eye and stimuli positions
            assert(length(EL_lag.s1_x{t})==length(EL_lag.x{b}{t}), 'stim/eye lengths don''t match!')
            assert(length(EL_lag.s1_y{t})==length(EL_lag.y{b}{t}), 'stim/eye lengths don''t match!')
        else
%             EL_lag.s2_x{t} = []; EL_lag.s2_y{t} = [];
            
            %interpolate the data to full length of uncut trial
            EL_lag.s2_x{t} = interp1(1:600, EL_lag.s2_x_600{t}, xi);
            EL_lag.s2_y{t} = interp1(1:600, EL_lag.s2_y_600{t}, xi);
            
            %remove the eye data (nan + blinks)
            EL_lag.s2_x{t}(EL_lag.to_remove{b}{t}) = []; EL_lag.s2_y{t}(EL_lag.to_remove{b}{t}) = [];
            
            %make sure we have the same length arrays for eye and stimuli positions
            assert(length(EL_lag.s2_x{t})==length(EL_lag.x{b}{t}), 'stim/eye lengths don''t match!')
            assert(length(EL_lag.s2_y{t})==length(EL_lag.y{b}{t}), 'stim/eye lengths don''t match!')
        end
        
    end
end


%% do the lag correlation to fix optimal lag

for b = 1:s   
for t = 1:Lag.tperS
    
    for n=1:500
        
        s1x_temp = EL_lag.s1_x{t}(1:end+1-n);
        s1y_temp = EL_lag.s1_y{t}(1:end+1-n);
        
        s2x_temp = EL_lag.s2_x{t}(1:end+1-n);
        s2y_temp = EL_lag.s2_y{t}(1:end+1-n);

        ex_temp = EL_lag.x{b}{t}(n:end);
        ey_temp = EL_lag.y{b}{t}(n:end);
        
        if b==1
            EL_lag.r_score_s1{t}(n,1) = corr2(s1x_temp, ex_temp);
            EL_lag.r_score_s1{t}(n,2) = corr2(s1y_temp, ey_temp);
        else
            EL_lag.r_score_s2{t}(n,1) = corr2(s2x_temp, ex_temp);
            EL_lag.r_score_s2{t}(n,2) = corr2(s2y_temp, ey_temp);
        end
    end
    
    if b==1
        EL_lag.max_corr_s1(t,1) = max(EL_lag.r_score_s1{t}(:,1)); % max correlation for x values
        EL_lag.max_corr_s1(t,2) = max(EL_lag.r_score_s1{t}(:,2)); % max correlation for y values
        
        EL_lag.max_lag_s1(t,1) = find(EL_lag.r_score_s1{t}(:,1)==max(EL_lag.r_score_s1{t}(:,1))); % optimal lag for x
        EL_lag.max_lag_s1(t,2) = find(EL_lag.r_score_s1{t}(:,2)==max(EL_lag.r_score_s1{t}(:,2))); % optimal lag for y
    else
        EL_lag.max_corr_s2(t,1) = max(EL_lag.r_score_s2{t}(:,1)); % max correlation for x values
        EL_lag.max_corr_s2(t,2) = max(EL_lag.r_score_s2{t}(:,2)); % max correlation for y values
        
        EL_lag.max_lag_s2(t,1) = find(EL_lag.r_score_s2{t}(:,1)==max(EL_lag.r_score_s2{t}(:,1))); % optimal lag for x
        EL_lag.max_lag_s2(t,2) = find(EL_lag.r_score_s2{t}(:,2)==max(EL_lag.r_score_s2{t}(:,2))); % optimal lag for y
    end
end
end


EL_lag.avgmaxlag.s1_x = round(mean(EL_lag.max_lag_s1(:,1)));
EL_lag.avgmaxlag.s1_y = round(mean(EL_lag.max_lag_s1(:,2)));
EL_lag.avgmaxlag.s2_x = round(mean(EL_lag.max_lag_s2(:,1)));
EL_lag.avgmaxlag.s2_y = round(mean(EL_lag.max_lag_s2(:,2)));

EL_lag.avgmaxlag.s1 = round(mean(mean(EL_lag.max_lag_s1())));
EL_lag.avgmaxlag.s2 = round(mean(mean(EL_lag.max_lag_s2())));
EL_lag.avgmaxlag.s1s2 = round(mean([EL_lag.avgmaxlag.s1, EL_lag.avgmaxlag.s2]));


disp('done!')
%% get the main messages

msgs_all = edf_main.Events.Messages.info;
msg_times_all = edf_main.Events.Messages.time;

%what messages did we pass on to Eyelink in the first place?
% sprintf('EVENT_Block%dTrial%d_Start', iBlock, iTrial)
% sprintf('EVENT_Block%dTrial%d_Fix', iBlock, iTrial)
% Eyelink('message', 'EVENT_TaskFaceOn');
% sprintf('EVENT_Block%dTrial%d_End', iBlock, iTrial)

e = 1;
for m = 1:length(msgs_all)
    if contains(msgs_all{m}, 'EVENT') %gather all the event messages (everything except calibration)
        evt_msgs{e} = msgs_all{m};
        evt_msgs_tm{e} = msg_times_all(m); %timestamps of each message
        e = e+1;
    end
end

e = 1;
for m = 1:length(evt_msgs)
    if contains(evt_msgs{m}, 'Task') %all messages relating to task
        task_msgs{e} = evt_msgs{m};
        task_msgs_tm{e} = evt_msgs_tm(m);
        e = e+1;
    end
end

%this will take some time to process
for b = Params.info.FroBlo:Params.info.To_Blo
    trial_msgs_tm{b} = {};
    for t = 1:Params.trialNum
        trial_msgs{b}{t} = {};
        trial_msgs_tm{b}{t} = [];
        for m = 1:length(evt_msgs)
            if contains(evt_msgs{m}, sprintf('_Block%dTrial%d_', b, t))
                trial_msgs{b}{t}{end+1} = evt_msgs{m};
                trial_msgs_tm{b}{t}(end+1) = evt_msgs_tm{m};
            end
        end
    end
end

messages.msgs_all = msgs_all;
messages.msg_times_all = msg_times_all;
messages.evt_msgs = evt_msgs;
messages.evt_msgs_tm = evt_msgs_tm;
messages.task_msgs = task_msgs;
messages.task_msgs_tm = task_msgs_tm;
messages.trial_msgs = trial_msgs;
messages.trial_msgs_tm = trial_msgs_tm;

% return msgs_all, msg_times_all, evt_msgs, evt_msgs_tm, task_msgs, task_msgs_tm, trial_msgs, trial_msgs_tm
%% calculate each trial separately

%collect the event messages from start to finish
%get block sequence of each subject (var = block)

%% call the functions

% for simplicity purposes: s1 = face, s2 = flower
face_EL = struct();
flower_EL = struct();
both_EL = struct();

face_EL = org_ELdata(face_EL, find(block==1), edf_main, Block, trial_msgs, trial_msgs_tm);
flower_EL = org_ELdata(flower_EL, find(block==2), edf_main, Block, trial_msgs, trial_msgs_tm);
both_EL = org_ELdata(both_EL, find(block==0), edf_main, Block, trial_msgs, trial_msgs_tm);
%%
face_EL = org_stimdata(face_EL, find(block==1), pT, pT2, Block);
flower_EL = org_stimdata(flower_EL, find(block==2), pT, pT2, Block);
both_EL = org_stimdata(both_EL, find(block==0), pT, pT2, Block);

%%
% face_EL = slidwin_ELcorr(face_EL, 1, EL_lag);
% flower_EL = slidwin_ELcorr(flower_EL, 1, EL_lag);
% both_EL = slidwin_ELcorr(both_EL, 1, EL_lag);

face_EL = slidwin_ELcorr(face_EL, 1, lag);
flower_EL = slidwin_ELcorr(flower_EL, 1, lag);
both_EL = slidwin_ELcorr(both_EL, 1, lag);

%%
face_EL = s1s2_corr_ext(face_EL, 1);
flower_EL = s1s2_corr_ext(flower_EL, 1);
both_EL = s1s2_corr_ext(both_EL, 1);

%%
function blocktype = org_ELdata(blocktype, whichblocks, edf_main, Block, trial_msgs, trial_msgs_tm)
% 1000 Hz and 60 Hz data are slightly offset - don't know why this discrepancy 
% exists BUT we are selecting the start/end points of 60Hz b/c it has 'lag' unlike
% with 1000 Hz = getting more data at the end by sacrificing the beginning

% whichblocks = find(block==att_id); %refs where each att-cond falls in expt

for bl = 1:length(whichblocks) % relative to each condition
    
    b = whichblocks(bl); % relative to exp sequence
    for t = 1:20
        
        %find where fixation ends and trial begins: 1000 Hz
        for m = 2:length(trial_msgs{b}{t})
            if contains(trial_msgs{b}{t}{m}, 'Fix')
%                 disp('still fixating...')
            else
%                 disp('stopped fixating!')
%                 disp(num2str(m))
                stopfix_id = m;
                break
            end
        end
        
        %find where fixation ends and trial begins: 60 Hz
        afterFix_idx = find(Block.EyeData{b}.FixDoneT{t} >= Block.flipTimes{b}{t}(2)); %which frames are after fixation
        endOfFix_id = afterFix_idx(1); %the first frame after fixation
        
        % find the start/end timestamp idx in Samples timeline based on 1000 Hz
        EL_trial.startmsg_abs_id{bl}{t} = find(edf_main.Samples.time==trial_msgs_tm{b}{t}(1)); %trial start before fixation
        EL_trial.startmsg_id{bl}{t} = find(edf_main.Samples.time==trial_msgs_tm{b}{t}(stopfix_id)); %trial begins after fixation
        EL_trial.endmsg_id{bl}{t}= find(edf_main.Samples.time==trial_msgs_tm{b}{t}(end)); %trial end
        
        % find the start/end timestamp idx in Samples timeline based on  60 Hz
        EL_trial.startTrial_abs_id{bl}{t} = find(edf_main.Samples.time==Block.EyeData{b}.time{t}(1)); %trial start before fixation
        EL_trial.startTrial_id{bl}{t} = find(edf_main.Samples.time==Block.EyeData{b}.time{t}(endOfFix_id)); %first point after fixation
        EL_trial.endTrial_id{bl}{t} = find(edf_main.Samples.time==Block.EyeData{b}.time{t}(end)); %end of trial
        
        % EL_trial.startmsg_id{bl}{t}-EL_trial.startTrial_id{bl}{t} %the discrepancy in starting points
        
        %the range of entire trial we Sample
        EL_trial.rangeTrial_id{bl}{t} = EL_trial.startTrial_id{bl}{t}:EL_trial.endTrial_id{bl}{t};
        EL_trial.rangeTrial_t{bl}{t} = edf_main.Samples.time(EL_trial.rangeTrial_id{bl}{t});
        EL_trial.trialDur_nocuts{bl}{t}= length(EL_trial.rangeTrial_id{bl}{t}); % in milliseconds, also the number of eye position points we start off with (not counting blinks)
        
        el.MISSING_DATA=-32768;
        % get all Eyelink values that are not missing data
        gx = edf_main.Samples.gx(edf_main.Samples.gx~=el.MISSING_DATA)';
        gy = edf_main.Samples.gy(edf_main.Samples.gy~=el.MISSING_DATA)';
        pa = edf_main.Samples.pa(edf_main.Samples.gy~=el.MISSING_DATA)';
        
        %get the values: keep a copy of the uncut values
        EL_x_full{bl}{t} = gx(EL_trial.rangeTrial_id{bl}{t});
        EL_y_full{bl}{t} = gy(EL_trial.rangeTrial_id{bl}{t});
        EL_a_full{bl}{t}= pa(EL_trial.rangeTrial_id{bl}{t});
        EL_t_full{bl}{t} = edf_main.Samples.time(EL_trial.rangeTrial_id{bl}{t})';
        
        %make a copy of values to edit so we have 2 versions: edited (below) and unedited (above)
        EL_x{bl}{t} = EL_x_full{bl}{t};
        EL_y{bl}{t} = EL_y_full{bl}{t};
        EL_a{bl}{t} = EL_a_full{bl}{t};
        EL_t{bl}{t} = EL_t_full{bl}{t};
        
        
        %% find the blinks
        
        % use Eyelink's built-in software to detect blinks, will remove the same
        % regions from x and y arrays
        EL_trial.blink_id{bl}{t} = find(edf_main.Events.Eblink.start>EL_trial.rangeTrial_t{bl}{t}(1) ...
            & edf_main.Events.Eblink.start<EL_trial.rangeTrial_t{bl}{t}(end)); %is there a blink that starts during trial
        %the 'corrupted' eye data before and after a blink to remove (ms)
        blink_buffer = 100;
        EL_trial.to_remove{bl}{t} = []; %add these indices to the ones we need to remove from the stimulus
        
        if ~isempty(EL_trial.blink_id{bl}{t})
            for blink = EL_trial.blink_id{bl}{t}
                blink_start_id = find(EL_trial.rangeTrial_t{bl}{t}==edf_main.Events.Eblink.start(blink)); % index of blink start time
                blink_start_t = edf_main.Events.Eblink.start(blink); % blink start time
                blink_end_id = find(EL_trial.rangeTrial_t{bl}{t}==edf_main.Events.Eblink.end(blink)); %index of blink end time
                blink_end_t = edf_main.Events.Eblink.end(blink); % blink end time
                
                % to make sure the buffer doesn't go beyond the trial bounds
                if blink_start_id-blink_buffer > 0
                    st_rem = blink_start_id-blink_buffer;
                else
                    st_rem = 1;
                end
                if blink_end_id+blink_buffer <= length(EL_t{bl}{t})
                    end_rem = blink_end_id+blink_buffer;
                else
                    end_rem = length(EL_t{bl}{t}); % x and y data will be same length
                end
                EL_trial.to_remove{bl}{t} = [EL_trial.to_remove{bl}{t} st_rem:end_rem]; %the 'corrupted' eye data before and after a blink to remove (ms)
                
                
            end
        end
        
        % to check that all x/y nans were included in removal
        x_isnan = isnan(EL_x_full{bl}{t});
        assert(sum(x_isnan==isnan(EL_y_full{bl}{t}))==length(x_isnan), 'NaNs don''t match!') %check that x/y nans match
        
        if sum(ismember(find(isnan(EL_x_full{bl}{t})), EL_trial.to_remove{bl}{t}))~=sum(isnan(EL_x_full{bl}{t}))
            %if perhaps the blink is in the beginning (or end) of the trial and Eyelink doesn't catch it?
            leftover_nan = setdiff(find(x_isnan), EL_trial.to_remove{bl}{t}); % x_isnan's NOT found in to_remove
            left_toremove = [];
            for l = leftover_nan
                left_toremove = [left_toremove l-blink_buffer:l+blink_buffer]; %add a buffer to each point
            end
            left_toremove = unique(left_toremove); %take unique points only (only the extremes will be kept)
            left_toremove = left_toremove(left_toremove>0 & left_toremove<length(EL_t{bl}{t})); %keep within bounds
            
            EL_trial.to_remove{bl}{t} = [EL_trial.to_remove{bl}{t} left_toremove]; %add to array to remove
            
        end
        
        assert(sum(ismember(find(x_isnan), EL_trial.to_remove{bl}{t}))==sum(x_isnan), sprintf('NaNs not deleted! block %d trial %d!', bl, t))
        
        
        % deletion: do it only once per trial or displacement will be wrong
        EL_x{bl}{t}(EL_trial.to_remove{bl}{t}) = [];
        EL_y{bl}{t}(EL_trial.to_remove{bl}{t}) = [];
        EL_a{bl}{t}(EL_trial.to_remove{bl}{t}) = [];
        EL_t{bl}{t}(EL_trial.to_remove{bl}{t}) = [];
        
        
       
    end %for-loop per trial       
end %for-loop per block

blocktype.EL_full.x = EL_x_full;
blocktype.EL_full.y = EL_y_full;
blocktype.EL_full.a = EL_a_full;
blocktype.EL_full.t = EL_t_full;

blocktype.EL_x = EL_x;
blocktype.EL_y = EL_y;
blocktype.EL_a = EL_a;
blocktype.EL_t = EL_t;

blocktype.EL_trial = EL_trial;

disp('all cleaned up!')
end % org_ELdata function call

%return: EL_x_full, EL_y_full, EL_a_full, EL_t_full, EL_x, EL_y, EL_a, EL_t, EL_trial(STRUCT)


%     d = zeros([1,length(EL_x_full{bl}{t})]);
%         d(EL_trial.to_remove{bl}{t}) = 600;
%         % figure
%         plot(EL_x{bl}{t})
%         hold on
%         plot(EL_x_full{bl}{t})
%         hold on
%         plot(d)
%         figure
%         plot(EL_y{bl}{t})
%         hold on
%         plot(EL_y_full{bl}{t})
%%

function blocktype = org_stimdata(blocktype, whichblocks, pT, pT2, Block)
% get the stimulus locations
% interpolate before removing invalid intervals (i.e. blinks)

%  whichblocks = find(block==att_id); %specific to each attend condition (out of 2)

for b = 1:length(whichblocks)
    bl = whichblocks(b);
    % original 4-coordinate stimuli positions across all trials in each block
    for iTrial = 1:20
        s1_pos{b}{iTrial} = pT{Block.t2Run{bl}(iTrial)}; %face
        s2_pos{b}{iTrial} = pT2{Block.t2Run{bl}(iTrial)}; %flower
    end  
end

% calculate average stimulus x and y values (center points) - (s1:face, s2:flower)
% ONLY gives us 600 points
for b = 1:length(whichblocks)
    s1_x_600{b} = {}; s1_y_600{b} = {};  
    s2_x_600{b} = {}; s2_y_600{b} = {};
    
    for trial = 1:20
        s1_x_600{b}{trial} = []; s1_y_600{b}{trial} = [];
        s2_x_600{b}{trial} = []; s2_y_600{b}{trial} = [];
        for i = 1:600 
            s1_x_600{b}{trial}(end+1)= mean([s1_pos{b}{trial}(i,1),s1_pos{(b)}{trial}(i,3)]);
            s1_y_600{b}{trial}(end+1)= mean([s1_pos{b}{trial}(i,2),s1_pos{(b)}{trial}(i,4)]);
            s2_x_600{b}{trial}(end+1)= mean([s2_pos{b}{trial}(i,1),s2_pos{(b)}{trial}(i,3)]);
            s2_y_600{b}{trial}(end+1)= mean([s2_pos{b}{trial}(i,2),s2_pos{(b)}{trial}(i,4)]);
        end
    end
end

% need to match the number of stimulus points to eye data points: linear interpolation
for b = 1:length(whichblocks)
    s1_x{b} = {}; s1_y{b} = {};  
    s2_x{b} = {}; s2_y{b} = {};
    
    for trial = 1:20
        % x and y should have the same number of points
        xi = linspace(1,600, blocktype.EL_trial.trialDur_nocuts{b}{trial}); % 1 to 'length(stim)', with 'length of uncut trial' amount of points
        
        s1_x{b}{trial} = interp1(1:600, s1_x_600{b}{trial}, xi);
        s1_y{b}{trial} = interp1(1:600, s1_y_600{b}{trial}, xi);
        s2_x{b}{trial} = interp1(1:600, s2_x_600{b}{trial}, xi);
        s2_y{b}{trial} = interp1(1:600, s2_y_600{b}{trial}, xi);
        
        %remove the eye data (nan + blinks)
        s1_x{b}{trial}(blocktype.EL_trial.to_remove{b}{trial}) = [];
        s1_y{b}{trial}(blocktype.EL_trial.to_remove{b}{trial}) = [];
        s2_x{b}{trial}(blocktype.EL_trial.to_remove{b}{trial}) = [];
        s2_y{b}{trial}(blocktype.EL_trial.to_remove{b}{trial}) = [];
    end
    
    %make sure we have the same length arrays for eye and stimuli positions
    assert(length(s1_x{b}{trial})==length(blocktype.EL_x{b}{trial}), 'stim/eye lengths don''t match!')
    assert(length(s1_y{b}{trial})==length(blocktype.EL_y{b}{trial}), 'stim/eye lengths don''t match!')
    assert(length(s2_x{b}{trial})==length(blocktype.EL_x{b}{trial}), 'stim/eye lengths don''t match!')
    assert(length(s2_y{b}{trial})==length(blocktype.EL_y{b}{trial}), 'stim/eye lengths don''t match!')
end

% return s1_pos, s2_pos, s1_x_600, s1_y_600, s2_x_600, s2_y_600, s1_x, s1_y, s2_x, s2_y

blocktype.s1_pos = s1_pos;
blocktype.s2_pos = s2_pos;

blocktype.s_600.s1x = s1_x_600;
blocktype.s_600.s1y = s1_y_600;
blocktype.s_600.s2x = s2_x_600;
blocktype.s_600.s2y = s2_y_600;

blocktype.s1_x = s1_x;
blocktype.s1_y = s1_y;
blocktype.s2_x = s2_x;
blocktype.s2_y = s2_y;

disp('stimuli are set!')
end % org_stimdata function call

%% sliding window of correlation: eye v stim

function blocktype = slidwin_ELcorr(blocktype, win_size, lag)
% EACH trial and EACH stim has an optimal lag
blocksPerCon = 2;
trialNo = 20;

win_frames = win_size*1000; %window size entered in seconds

for b = 1:blocksPerCon
    
    eye_x = blocktype.EL_x{b};
    eye_y = blocktype.EL_y{b};
    
    stim1_x = blocktype.s1_x{b};
    stim1_y = blocktype.s1_y{b};
    
    stim2_x = blocktype.s2_x{b};
    stim2_y = blocktype.s2_y{b};
    
    for t = 1:trialNo
        
        if isfield(lag.avgmaxlag, 'facehouse')
            lag.avgmaxlag.s1s2 = lag.avgmaxlag.facehouse;
        end

        %1-1 lineup of eye-x and stim-1 and stim-2
        x1_trial = eye_x{t}(lag.avgmaxlag.s1s2:end);
        x2_trial = eye_x{t}(lag.avgmaxlag.s1s2:end);
        s1x_trial = stim1_x{t}(1:end-lag.avgmaxlag.s1s2+1);
        s2x_trial = stim2_x{t}(1:end-lag.avgmaxlag.s1s2+1);

        y1_trial = eye_y{t}(lag.avgmaxlag.s1s2:end);
        y2_trial = eye_y{t}(lag.avgmaxlag.s1s2:end);
        s1y_trial = stim1_y{t}(1:end-lag.avgmaxlag.s1s2+1);
        s2y_trial = stim2_y{t}(1:end-lag.avgmaxlag.s1s2+1);

        win_shifts_s1x = length(s1x_trial)-win_frames+1;
        win_shifts_s2x = length(s2x_trial)-win_frames+1;
        win_shifts_s1y = length(s1y_trial)-win_frames+1;
        win_shifts_s2y = length(s2y_trial)-win_frames+1;

        corr_win.s1.x{t} = zeros([1, win_shifts_s1x]);
        corr_win.s1.y{t} = zeros([1, win_shifts_s1y]);
        corr_win.s2.x{t} = zeros([1, win_shifts_s2x]);
        corr_win.s2.y{t} = zeros([1, win_shifts_s2y]);

        for w = 1:win_shifts_s1x        
            corr_win.s1.x{t}(w) = corr2(x1_trial(w:w+win_frames-1), s1x_trial(w:w+win_frames-1));
        end
        for w = 1:win_shifts_s2x        
            corr_win.s2.x{t}(w) = corr2(x2_trial(w:w+win_frames-1), s2x_trial(w:w+win_frames-1));
        end
        for w = 1:win_shifts_s1y        
            corr_win.s1.y{t}(w) = corr2(y1_trial(w:w+win_frames-1), s1y_trial(w:w+win_frames-1));
        end
        for w = 1:win_shifts_s2y        
            corr_win.s2.y{t}(w) = corr2(y2_trial(w:w+win_frames-1), s2y_trial(w:w+win_frames-1));
        end

        corrAnalyses.eye_x1{b}{t} = x1_trial; %eye-x values shifted by optimal lag for stim1
        corrAnalyses.eye_x2{b}{t} = x2_trial; %eye-x values shifted by optimal lag for stim2

        corrAnalyses.eye_y1{b}{t} = y1_trial; %eye-y values shifted by optimal lag for stim1
        corrAnalyses.eye_y2{b}{t} = y2_trial; %eye-y values shifted by optimal lag for stim2

        corrAnalyses.s1x{b}{t} = s1x_trial;
        corrAnalyses.s1y{b}{t} = s1y_trial;

        corrAnalyses.s2x{b}{t} = s2x_trial;
        corrAnalyses.s2y{b}{t} = s2y_trial;

    end

    % the correlation results for each sliding window
    blocktype.corr_win{b} = corr_win; 
    
    % the mean correlation for each trial (1,:) = x, (2,:) = y
    blocktype.mean_s1{b} = zeros([2, trialNo]);
    blocktype.mean_s2{b} = zeros([2, trialNo]);
    for t = 1:trialNo
        blocktype.corr_mean.s1{b}(1,t) = mean(corr_win.s1.x{t}); 
        blocktype.corr_mean.s1{b}(2,t) = mean(corr_win.s1.y{t});
        blocktype.corr_mean.s2{b}(1,t) = mean(corr_win.s2.x{t});
        blocktype.corr_mean.s2{b}(2,t) = mean(corr_win.s2.y{t});
    end

end

%the inputs for each sliding window correlation output, shifted by the lag
blocktype.corr_inputs = corrAnalyses;

%calculate the average correlation for each stimulus across blocks
blocktype.corr_mean.s1x = mean([blocktype.corr_mean.s1{1}(1,:) blocktype.corr_mean.s1{2}(1,:)]);
blocktype.corr_mean.s1y = mean([blocktype.corr_mean.s1{1}(2,:) blocktype.corr_mean.s1{2}(2,:)]);
blocktype.corr_mean.s2x = mean([blocktype.corr_mean.s2{1}(1,:) blocktype.corr_mean.s2{2}(1,:)]);
blocktype.corr_mean.s2y = mean([blocktype.corr_mean.s2{1}(2,:) blocktype.corr_mean.s2{2}(2,:)]);

%calculate the average correlation for each stimulus across x/y    
blocktype.corr_mean_s1 = mean([blocktype.corr_mean.s1x, blocktype.corr_mean.s1y]);
blocktype.corr_mean_s2 = mean([blocktype.corr_mean.s2x, blocktype.corr_mean.s2y]);

disp(['avg face corr: ' num2str(blocktype.corr_mean_s1)])
disp(['avg flower corr: ' num2str(blocktype.corr_mean_s2)])
end

 %% temporal correlation: stim1(face) v stim 2(house)

function blocktype = s1s2_corr_ext(blocktype, win_size)
%calculate stim1-stim2 correlation using the interpolated position values = extended array

blocksPerCon = 2;
trialNo = 20;

win_frames = win_size *1000; %window size entered in seconds


for b = 1:blocksPerCon
    
    stim1_x = blocktype.s1_x{b};
    stim2_x = blocktype.s2_x{b};
    stim1_y = blocktype.s1_y{b};
    stim2_y = blocktype.s2_y{b};
    
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
    
    % the correlation results for each sliding window
    blocktype.corr_win{b}.s1s2_x = corr_s1s2x;
    blocktype.corr_win{b}.s1s2_y = corr_s1s2y;
    
    % average correlation results for each trial
    blocktype.corr_mean.s1s2{b} = zeros([2, trialNo]);
    for t = 1:trialNo
        blocktype.corr_mean.s1s2{b}(1,t) = mean(corr_s1s2x{t});
        blocktype.corr_mean.s1s2{b}(2,t) = mean(corr_s1s2y{t});
    end
  
end

% average correlation across blocks
blocktype.corr_mean.s1s2x = mean([blocktype.corr_mean.s1s2{1}(1,:) blocktype.corr_mean.s1s2{2}(1,:)]);
blocktype.corr_mean.s1s2y =  mean([blocktype.corr_mean.s1s2{1}(2,:) blocktype.corr_mean.s1s2{2}(2,:)]);

%average correlation across x and y
blocktype.corr_mean.s1s2xy = mean([blocktype.corr_mean.s1s2x blocktype.corr_mean.s1s2y]);   

disp(['avg face-flower corr: ' num2str(blocktype.corr_mean.s1s2xy)])
disp('done!')
end
