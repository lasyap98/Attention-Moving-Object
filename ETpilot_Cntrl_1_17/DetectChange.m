function [task_frames, change_num] = DetectChange(change_duration)

global Display Params
%% create perceptual task
% for one specific trial, change both

% rng('shuffle'); %can't run this on ET computer in 425
% Display.refresh = 1/60;
lagOn = 0;
%account for both stimulus changes to prevent overlap
if Params.info.Prac || lagOn
    change_num = 10; %want all practice trials to have 10/2 = 5 events for each stimulus
else
    change_num = 2*randsample(3:5, 1); %can have between 3-5 changes in a trial, prevent predictability 
end

Params.change_dur = change_duration; %sec
Params.change_frames = round(Params.change_dur/Display.refresh);
poss_chg_no = round(Params.trialLength/Params.change_dur); %possible number of change events in a 10sec trial

gap = round(.6/Params.change_dur); % sec
%% simpler method?

trial_evts = zeros([1,poss_chg_no]);
while (sum(trial_evts)~=change_num)
    trial_evts = zeros([1,poss_chg_no]);
    p = 2;
    while p <= poss_chg_no-1
        evt_success = rand();
        if (evt_success >.5) && (sum(trial_evts)<change_num)
            trial_evts(p) = 1;
            p = p+gap+1;
        else
            p = p+1;
        end
    end
end
% sum(trial_evts)

%%
%to prevent overlap of stimulus events
event_idx = find(trial_evts==1); % both stimulus change events
s1_evt_idx = randsample(event_idx, change_num/2); % randomly assign half of these events to one stimulus 
s2_evt_idx = setdiff(event_idx, s1_evt_idx); % assign the rest of these events to the other stimulus 

s1_trial_evts = trial_evts;
s1_trial_evts(s2_evt_idx) = 0; %events in s2 shouldn't overlap with s1

s2_trial_evts = trial_evts;
s2_trial_evts(s1_evt_idx) = 0; %events in s1 shouldn't overlap with s2

% expand the events by change_frames to span the length of a trial
% task_frames(1,:) = repelem(s1_trial_evts,Params.change_frames);
% task_frames(2,:) = repelem(s2_trial_evts,Params.change_frames);

task_frames_mat = zeros([length(s1_trial_evts),Params.change_frames]);
task_frames_mat2 = zeros([length(s2_trial_evts),Params.change_frames]);
tf1 = []; tf2 = [];
task_frames = [];
for t = 1:length(s2_trial_evts)
    for i = 1:Params.change_frames
        task_frames_mat(t,i) = s1_trial_evts(t);
        task_frames_mat2(t,i) = s2_trial_evts(t);
    end
    tf1 = [tf1 task_frames_mat(t,:)];
    tf2 = [tf2 task_frames_mat2(t,:)];
end

task_frames(1,:) = tf1;
task_frames(2,:) = tf2;



end


