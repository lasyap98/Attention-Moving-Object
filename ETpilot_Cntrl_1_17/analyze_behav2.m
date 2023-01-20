 % RT should only count if they fall within a certain time-frame after change
 % my average was around 500ms, set to 1 sec?
 resp_interval = [.1,1]; %sec
 resp_interval_frames = round(resp_interval/refresh); %frames
  
 if Params.sub.subject==1
     block = [1, 2, 0, 1, 2, 0]; %face, house
 elseif Params.sub.subject==2 
     block = [0, 1, 2, 0, 1, 2]; %both, face, house
 elseif Params.sub.subject==3
     block = [2, 0, 1, 2, 0, 1];
 end


face.eyedata = Block.EyeData(block==1);
face.attend = attend(block==1);
face.changes = Changes_block(block==1);
face.fliptimes = Block.flipTimes(block==1);
face.keyi = Block.key_i(block==1);
face.keyCode_trial = Block.keyCode_trial(block==1);
face.keySecs_trial = Block.keySecs_trial(block==1);
face.taskframes = Task_block(block==1);

%%
house.eyedata = Block.EyeData(block==2);
house.attend = attend(block==2);
house.changes = Changes_block(block==2);
house.fliptimes = Block.flipTimes(block==2);
house.keyi = Block.key_i(block==2);
house.keyCode_trial = Block.keyCode_trial(block==2);
house.keyCode_trial = Block.keyCode_trial(block==2);
house.keySecs_trial = Block.keySecs_trial(block==2);
house.taskframes = Task_block(block==2);

%%
both.attend = attend(block==0);
both.eyedata = Block.EyeData(block==0);
both.changes = Changes_block(block==0);
both.fliptimes = Block.flipTimes(block==0);
both.keyi = Block.key_i(block==0);
both.keyCode_trial = Block.keyCode_trial(block==0);
both.keySecs_trial = Block.keySecs_trial(block==0);
both.taskframes = Task_block(block==0);

%% start stop idx of task - house
for b = 1:length(house.taskframes)
    house_idxtask{b} = {};
    house_idxstart{b} = {};
    house_idxstop{b} = {};
    for t = 1:length(house.taskframes{b})
        house_idxtask{b}{t} = [];
        house_idxstart{b}{t} = [];
        house_idxstop{b}{t} = {};
        for f = 1:length(house.taskframes{b}{t})
            if house.taskframes{b}{t}(2,f)==1
                house_idxtask{b}{t}(end+1) = f; %all the frames the task was on
            end
        end
        len = length(house_idxtask{b}{t});
        house_idxstart{b}{t} = [house_idxtask{b}{t}(1), house_idxtask{b}{t}(Params.change_frames+1), house_idxtask{b}{t}((Params.change_frames*2)+1)];
        house_idxstop{b}{t} = [house_idxtask{b}{t}(Params.change_frames), house_idxtask{b}{t}(Params.change_frames*2), house_idxtask{b}{t}(Params.change_frames*3)];
        if (len/Params.change_frames)==4
            house_idxstart{b}{t}(end+1) = house_idxtask{b}{t}((Params.change_frames*3)+1);
            house_idxstop{b}{t}(end+1) = house_idxtask{b}{t}(Params.change_frames*4);
        elseif (len/Params.change_frames)==5
            house_idxstart{b}{t}(end+1) = house_idxtask{b}{t}((Params.change_frames*3)+1);
            house_idxstart{b}{t}(end+1) = house_idxtask{b}{t}((Params.change_frames*4)+1);
            
            house_idxstop{b}{t}(end+1) = house_idxtask{b}{t}(Params.change_frames*4);
            house_idxstop{b}{t}(end+1) = house_idxtask{b}{t}(Params.change_frames*5);
        end

    end
end

house.task.startIdx = house_idxstart;
house.task.stopIdx = house_idxstop;
house.task.taskIdx = house_idxtask;
%% calculate valid intervals for correct response
% by frame
for b = 1:length(house_idxstop)
    house.poss_resp_i{b} = {};
    for t = 1:length(house_idxstop{b})
        house.poss_resp_i{b}{t} = [];
        for i = 1:length(house_idxstop{b}{t})
            %can make resp starting from 100ms after stim starts, up til 1000ms after stim ends 
            house.poss_resp_i{b}{t} = [house.poss_resp_i{b}{t} house_idxstart{b}{t}(i)+resp_interval_frames(1):house_idxstop{b}{t}(i) + resp_interval_frames(2)];
        end
        house.poss_resp_i{b}{t} = sort(house.poss_resp_i{b}{t});
    end
end

%by time

for b = 1:length(house_idxstop)
    for t = 1:length(house_idxstop{b})
        house.task.stopTimes{b}{t} = house.fliptimes{b}{t}(house_idxstop{b}{t}+1); %end after it turns off, which is one extra frame after
        house.task.startTimes{b}{t} = house.fliptimes{b}{t}(house_idxstart{b}{t});
        for i = 1:length(house.task.startTimes{b}{t})
            house.possRT{b}{t}(i,:) = [house.task.startTimes{b}{t}(i)+.1 house.task.stopTimes{b}{t}(i)+1];
        end
    end
end

%% accuracy - house
house.correct = 0;
house.correct_sing = 0;
house.correct_doub = 0;
house.doub_count = 0;
house.sing_count = 0;
house.total_count = 0;
house.misses = 0;

for b = 1:length(house_idxstop)
    for t = 1:length(house_idxstop{b})
        
        trial_keys = house.keyi{b}{t};
        trial_keySecs = house.keySecs_trial{b}{t}(house.keySecs_trial{b}{t}~=0);
        
        %initialize trial counts
        house.trial_correct{b}{t} = 0; %correct specific to the trial
        house.trial_falsealarms{b}{t} = 0; %false alarms specific to the trial
        house.trial_miss{b}{t} = 0; %misses specific to the trial
        house.correct_keys{b}{t} = []; %which keys were correct specific to trial
        house.RT.secs{b}{t} = []; %RT in actual seconds specific to trial
        house.RT.frames{b}{t} = []; %RT in frames specific to trial
        house.resp_summary{b}{t} = zeros([1,house.changes{b}{t}/2]); %a lineup of each stim event, to be noted whether correct or not
        
        for i = 1:length(house_idxstop{b}{t})
            % don't count the first .1 secs after stimulus starts 
            evt_resp_interval = house_idxstart{b}{t}(i)+resp_interval_frames(1):house_idxstop{b}{t}(i)+resp_interval_frames(2);
            [probe_resp, resp_id] = intersect(trial_keys, evt_resp_interval);
            
            if ~isempty(resp_id) %if there is a hit
                
                %select only the first valid response
                RT = trial_keySecs(resp_id(1))-house.task.startTimes{b}{t}(i); %RT in actual seconds
                RT_id = probe_resp(1)-house_idxstart{b}{t}(i); %RT in frames 

                house.trial_correct{b}{t} = house.trial_correct{b}{t} +1; %correct specific to the trial
                house.correct_keys{b}{t} = [house.correct_keys{b}{t} probe_resp(1)]; %note down frame that correct key was pressed
                house.RT.secs{b}{t} = [house.RT.secs{b}{t} RT];
                house.RT.frames{b}{t} = [house.RT.frames{b}{t} RT_id];
                house.resp_summary{b}{t}(i) = 1; %note that it was correct in the lineup
                
                %remove that key from list so we don't count it twice
                trial_keys(resp_id(1)) = []; 
                trial_keySecs(resp_id(1)) = []; 

            else %it's a miss
                house.trial_miss{b}{t} = house.trial_miss{b}{t} + 1;
            end
        end
        

        %if there are any leftover keys from the trial that didn't count towards a hit, then it's a false alarm
        house.trial_falsealarms{b}{t} = length(trial_keys);

        %calculate FA-ratio
        house.FA_ratio{b}{t} = house.trial_falsealarms{b}{t}/length(house.keyi{b}{t}); % #false alarms/total responses
               
        %update count
        if house.attend{b}{t}==2 %on trials where only the attended stimulus changes
            house.correct_sing = house.correct_sing + house.trial_correct{b}{t};
            house.sing_count = house.sing_count + house.changes{b}{t}/2;
        elseif house.attend{b}{t}==3 %on trials where both the stimuli changes
            house.correct_doub = house.correct_doub + house.trial_correct{b}{t};
            house.doub_count = house.doub_count + house.changes{b}{t}/2;
        end
        
        % consolidate counts: correct + misses = total
        house.correct = house.correct + house.trial_correct{b}{t}; %add each trial's correct responses to total count
        house.misses = house.misses + house.trial_miss{b}{t};
        house.total_count = house.total_count+ house.changes{b}{t}/2;
        
        
        
    end
end
assert(house.correct+house.misses==house.total_count, 'house counts don''t match up!')

house.accuracy = house.correct/house.total_count;
house.mean_FAR = nanmean([house.FA_ratio{1}{:} house.FA_ratio{2}{:}]);
disp(['house accuracy: ', num2str(house.accuracy)]);      
disp(['house FAR: ', num2str(house.mean_FAR)]);      
%% start, stop, idx of task for face
for b = 1:length(face.taskframes)
    face_idxtask{b} = {};
    face_idxstart{b} = {};
    face_idxstop{b} = {};
    for t = 1:length(face.taskframes{b})
        face_idxtask{b}{t} = [];
        face_idxstart{b}{t} = [];
        face_idxstop{b}{t} = {};
        for f = 1:length(face.taskframes{b}{t})
            if face.taskframes{b}{t}(1,f)==1
                face_idxtask{b}{t}(end+1) = f;
            end
        end
        len = length(face_idxtask{b}{t});
        face_idxstart{b}{t} = [face_idxtask{b}{t}(1), face_idxtask{b}{t}(Params.change_frames+1), face_idxtask{b}{t}((Params.change_frames*2)+1)];
        face_idxstop{b}{t} = [face_idxtask{b}{t}(Params.change_frames), face_idxtask{b}{t}(Params.change_frames*2), face_idxtask{b}{t}(Params.change_frames*3)];
        if (len/Params.change_frames)==4
            face_idxstart{b}{t}(end+1) = face_idxtask{b}{t}((Params.change_frames*3)+1);
            face_idxstop{b}{t}(end+1) = face_idxtask{b}{t}(Params.change_frames*4);
        elseif (len/Params.change_frames)==5
            face_idxstart{b}{t}(end+1) = face_idxtask{b}{t}((Params.change_frames*3)+1);
            face_idxstart{b}{t}(end+1) = face_idxtask{b}{t}((Params.change_frames*4)+1);
            
            face_idxstop{b}{t}(end+1) = face_idxtask{b}{t}(Params.change_frames*4);
            face_idxstop{b}{t}(end+1) = face_idxtask{b}{t}(Params.change_frames*5);
        end
    end
end

face.task.startIdx = face_idxstart;
face.task.stopIdx = face_idxstop;
face.task.taskIdx = face_idxtask;
%% calculate valid intervals for correct response
%by frames
for b = 1:length(face_idxstop)
    face.poss_resp_i{b} = {};
    for t = 1:length(face_idxstop{b})
        face.poss_resp_i{b}{t} = [];
        for i = 1:length(face_idxstop{b}{t})
            %can make resp starting from 100ms after stim starts, up til 1000ms after stim ends
            face.poss_resp_i{b}{t} =  [face.poss_resp_i{b}{t} face_idxstart{b}{t}(i)+resp_interval_frames(1) : face_idxstop{b}{t}(i) + resp_interval_frames(2)];
        end
        face.poss_resp_i{b}{t} = sort(face.poss_resp_i{b}{t});
    end
end

%by time

for b = 1:length(face_idxstop)
    for t = 1:length(face_idxstop{b})
        face.task.stopTimes{b}{t} = face.fliptimes{b}{t}(face_idxstop{b}{t}+1); %end after it turns off, which is one extra frame after
        face.task.startTimes{b}{t} = face.fliptimes{b}{t}(face_idxstart{b}{t});
        for i = 1:length(face.task.startTimes{b}{t})
            face.possRT{b}{t}(i,:) = [face.task.startTimes{b}{t}(i)+.1 face.task.stopTimes{b}{t}(i)+1];
        end
    end
end
%% calculate accuracy - face
face.correct = 0;
face.correct_sing = 0;
face.correct_doub = 0;
face.doub_count = 0;
face.sing_count = 0;
face.total_count = 0;
face.misses = 0;

for b = 1:length(face_idxstop)
    for t = 1:length(face_idxstop{b})
        
        trial_keys = face.keyi{b}{t};
        trial_keySecs = face.keySecs_trial{b}{t}(face.keySecs_trial{b}{t}~=0);
        
        assert(length(trial_keys)==length(trial_keySecs), 'key_i not the same length as key_Secs!!')
        
        %initialize trial counts
        face.trial_correct{b}{t} = 0; %correct specific to the trial
        face.trial_falsealarms{b}{t} = 0; %false alarms specific to the trial
        face.trial_miss{b}{t} = 0; %misses specific to the trial
        face.correct_keys{b}{t} = []; %which keys were correct specific to trial
        face.RT.secs{b}{t} = []; %RT in actual seconds specific to trial
        face.RT.frames{b}{t} = []; %RT in frames specific to trial
        face.resp_summary{b}{t} = zeros([1,face.changes{b}{t}/2]); %a lineup of each stim event, to be noted whether correct or not
        
        for i = 1:length(face_idxstop{b}{t})
            % don't count the first .1 secs after stimulus starts 
            evt_resp_interval = face_idxstart{b}{t}(i)+resp_interval_frames(1):face_idxstop{b}{t}(i)+resp_interval_frames(2);
            [probe_resp, resp_id] = intersect(trial_keys, evt_resp_interval);
            
            if ~isempty(resp_id) %if there is a hit
                
                %take only the first valid response
                RT = trial_keySecs(resp_id(1))-face.task.startTimes{b}{t}(i); %RT in actual seconds
                RT_id = probe_resp(1)-face_idxstart{b}{t}(i); %RT in frames 

                face.trial_correct{b}{t} = face.trial_correct{b}{t} +1; %correct specific to the trial
                face.correct_keys{b}{t} = [face.correct_keys{b}{t} probe_resp(1)]; %note down frame that correct key was pressed
                face.RT.secs{b}{t} = [face.RT.secs{b}{t} RT];
                face.RT.frames{b}{t} = [face.RT.frames{b}{t} RT_id];
                face.resp_summary{b}{t}(i) = 1; %note that correct in lineup
                
                %remove that key from list so we don't count it twice
                trial_keys(resp_id(1)) = []; 
                trial_keySecs(resp_id(1)) = []; 

            else %it's a miss
                face.trial_miss{b}{t} = face.trial_miss{b}{t} + 1;
            end
        end
        
        %if there are any leftover keys from the trial that didn't count towards a hit, then it's a false alarm
        face.trial_falsealarms{b}{t} = length(trial_keys);

        %calculate FA-ratio
        face.FA_ratio{b}{t} = face.trial_falsealarms{b}{t}/length(face.keyi{b}{t}); % #false alarms/total responses
               
        %update count
        if face.attend{b}{t}==1 %on trials where only the attended stimulus changes
            face.correct_sing = face.correct_sing + face.trial_correct{b}{t};
            face.sing_count = face.sing_count + face.changes{b}{t}/2;
        elseif face.attend{b}{t}==3 %on trials where both the stimuli changes
            face.correct_doub = face.correct_doub + face.trial_correct{b}{t};
            face.doub_count = face.doub_count + face.changes{b}{t}/2;
        end
        
        % consolidate counts: correct + misses = total
        face.correct = face.correct + face.trial_correct{b}{t}; %add each trial's correct responses to total count
        face.misses = face.misses + face.trial_miss{b}{t};
        face.total_count = face.total_count+ face.changes{b}{t}/2;
        
        
        
    end
end
assert(face.correct+face.misses==face.total_count, 'face counts don''t match up!')

face.accuracy = face.correct/face.total_count;
face.mean_FAR = nanmean([face.FA_ratio{1}{:} face.FA_ratio{2}{:}]);
disp(['face accuracy: ', num2str(face.accuracy)]);      
disp(['face FAR: ', num2str(face.mean_FAR)]);    
%% start stop idx of task - both
for b = 1:length(both.taskframes)
   
    both_faceidxstart{b} = {};
    both_faceidxstop{b} = {};
    
    both_houseidxstart{b} = {};
    both_houseidxstop{b} = {};
    
    for t = 1:length(both.taskframes{b})
       
        both_faceidxstart{b}{t} = [];
        both_faceidxstop{b}{t} = [];
        
        both_houseidxstart{b}{t} = [];
        both_houseidxstop{b}{t} = [];
        
        both_faceidxtask{b}{t} = find(both.taskframes{b}{t}(1,:)==1);
        both_houseidxtask{b}{t} = find(both.taskframes{b}{t}(2,:)==1);
        
        len = length(both_houseidxtask{b}{t});
        both_houseidxstart{b}{t} = [both_houseidxtask{b}{t}(1), both_houseidxtask{b}{t}(Params.change_frames+1), both_houseidxtask{b}{t}((Params.change_frames*2)+1)];
        both_houseidxstop{b}{t} = [both_houseidxtask{b}{t}(Params.change_frames), both_houseidxtask{b}{t}(Params.change_frames*2), both_houseidxtask{b}{t}(Params.change_frames*3)];
        
        
        if (len/Params.change_frames)==4
            both_houseidxstart{b}{t}(end+1) = both_houseidxtask{b}{t}((Params.change_frames*3)+1);
            both_houseidxstop{b}{t}(end+1) = both_houseidxtask{b}{t}(Params.change_frames*4);
        elseif (len/Params.change_frames)==5
            both_houseidxstart{b}{t}(end+1) = both_houseidxtask{b}{t}((Params.change_frames*3)+1);
            both_houseidxstart{b}{t}(end+1) = both_houseidxtask{b}{t}((Params.change_frames*4)+1);
            
            both_houseidxstop{b}{t}(end+1) = both_houseidxtask{b}{t}(Params.change_frames*4);
            both_houseidxstop{b}{t}(end+1) = both_houseidxtask{b}{t}(Params.change_frames*5);
        end
        
        len = length(both_faceidxtask{b}{t});
        both_faceidxstart{b}{t} = [both_faceidxtask{b}{t}(1), both_faceidxtask{b}{t}(Params.change_frames+1), both_faceidxtask{b}{t}((Params.change_frames*2)+1)];
        both_faceidxstop{b}{t} = [both_faceidxtask{b}{t}(Params.change_frames), both_faceidxtask{b}{t}(Params.change_frames*2), both_faceidxtask{b}{t}(Params.change_frames*3)];
        
        if (len/Params.change_frames)==4
            both_faceidxstart{b}{t}(end+1) = both_faceidxtask{b}{t}((Params.change_frames*3)+1);
            both_faceidxstop{b}{t}(end+1) = both_faceidxtask{b}{t}(Params.change_frames*4);
        elseif (len/Params.change_frames)==5
            both_faceidxstart{b}{t}(end+1) = both_faceidxtask{b}{t}((Params.change_frames*3)+1);
            both_faceidxstart{b}{t}(end+1) = both_faceidxtask{b}{t}((Params.change_frames*4)+1);
            
            both_faceidxstop{b}{t}(end+1) = both_faceidxtask{b}{t}(Params.change_frames*4);
            both_faceidxstop{b}{t}(end+1) = both_faceidxtask{b}{t}(Params.change_frames*5);
        end
        
    end
end


both.task.house_startIdx = both_houseidxstart;
both.task.house_stopIdx = both_houseidxstop;
both.task.house_taskIdx = both_houseidxtask;

both.task.face_startIdx = both_faceidxstart;
both.task.face_stopIdx = both_faceidxstop;
both.task.face_taskIdx = both_faceidxtask;

%% calculate valid intervals for correct response - both

%by frames
for b = 1:length(both_houseidxstop)
    both.houseposs_resp_i{b} = {};
    for t = 1:length(both_houseidxstop{b})
        both.houseposs_resp_i{b}{t} = both_houseidxtask{b}{t};
        for i = 1:length(both_houseidxstop{b}{t})
            both.houseposs_resp_i{b}{t} = [both.houseposs_resp_i{b}{t} both_houseidxstop{b}{t}(i)+1:both_houseidxstop{b}{t}(i) + resp_interval_frames(2)];
        end
        both.houseposs_resp_i{b}{t} = sort(both.houseposs_resp_i{b}{t});
    end
end

for b = 1:length(both_faceidxstop)
    both.faceposs_resp_i{b} = {};
    for t = 1:length(both_faceidxstop{b})
        both.faceposs_resp_i{b}{t} = both_faceidxtask{b}{t};
        for i = 1:length(both_faceidxstop{b}{t})
            both.faceposs_resp_i{b}{t} = [both.faceposs_resp_i{b}{t} both_faceidxstop{b}{t}(i)+1:both_faceidxstop{b}{t}(i) + resp_interval_frames(2)];
        end
        both.faceposs_resp_i{b}{t} = sort(both.faceposs_resp_i{b}{t});
    end
end


%by time
for b = 1:length(both_houseidxstop)
    for t = 1:length(both_houseidxstop{b})
        both.task.house_stopTimes{b}{t} = both.fliptimes{b}{t}(both_houseidxstop{b}{t}+1); %end after it turns off, which is one extra frame after
        both.task.house_startTimes{b}{t} = both.fliptimes{b}{t}(both_houseidxstart{b}{t});
        for i = 1:length(both.task.house_startTimes{b}{t})
            both.housepossRT{b}{t}(i,:) = [both.task.house_startTimes{b}{t}(i)+.1 both.task.house_stopTimes{b}{t}(i)+1];
        end
    end
end

for b = 1:length(both_faceidxstop)
    for t = 1:length(both_faceidxstop{b})
        both.task.face_stopTimes{b}{t} = both.fliptimes{b}{t}(both_faceidxstop{b}{t}+1); %end after it turns off, which is one extra frame after
        both.task.face_startTimes{b}{t} = both.fliptimes{b}{t}(both_faceidxstart{b}{t});
        for i = 1:length(both.task.face_startTimes{b}{t})
            both.facepossRT{b}{t}(i,:) = [both.task.face_startTimes{b}{t}(i)+.1 both.task.face_stopTimes{b}{t}(i)+1];
        end
    end
end

%% accuracy - both
both.facecorrect = 0;
both.housecorrect = 0;
both.face_count = 0;
both.house_count = 0;

for b = 1:length(both_faceidxstart)
    for t = 1:length(both_faceidxstart{b})
        [all_start_idx, start_sort] = sort([both_faceidxstart{b}{t}, both_houseidxstart{b}{t}]);
        
        all_stop_idx = [both_faceidxstop{b}{t}, both_houseidxstop{b}{t}];
        all_stop_idx = all_stop_idx(start_sort);
        
        all_start_times = [both.task.face_startTimes{b}{t}, both.task.house_startTimes{b}{t}];
        all_start_times = all_start_times(start_sort);
        
        %the order of events based on task
        task_start = [ones([1,length(both_faceidxstart{b}{t})]) 2*ones([1,length(both_houseidxstart{b}{t})])];
        task_start = task_start(start_sort); 
        
        trial_keys = both.keyi{b}{t};
        trial_keySecs = both.keySecs_trial{b}{t}(both.keySecs_trial{b}{t}~=0);
        
        %initialize trial counts
        both.trial_correct.face{b}{t} = 0;
        both.trial_correct.house{b}{t} = 0;
        both.trial_falsealarm{b}{t} = 0;
        both.correct_keys.face{b}{t} = []; 
        both.correct_keys.house{b}{t} = []; 
        both.face_RT.secs{b}{t} = [];
        both.face_RT.frames{b}{t} = [];
        both.house_RT.secs{b}{t} = [];
        both.house_RT.frames{b}{t} = [];
        both.face_resp_summary{b}{t} = zeros([1,both.changes{b}{t}/2]); %a lineup of each stim event, to be noted whether correct or not
        both.house_resp_summary{b}{t} = zeros([1,both.changes{b}{t}/2]); %a lineup of each stim event, to be noted whether correct or not
        
        for i = 1:length(all_start_idx)
            % don't count the first .1 secs after stimulus starts 
            evt_resp_interval = all_start_idx(i)+resp_interval_frames(1):all_stop_idx(i)+resp_interval_frames(2);
            [probe_resp, resp_id] = intersect(trial_keys, evt_resp_interval);
            
            if ~isempty(resp_id) %if there is a valid response
                
                %take only the first valid response
                RT = trial_keySecs(resp_id(1))-all_start_times(i); %RT in actual seconds
                RT_id = probe_resp(1)-all_start_idx(i); %RT in frames
                
                if task_start(i) == 1 %if this was a face task
                    
                    both.facecorrect = both.facecorrect + 1;
                    both.trial_correct.face{b}{t} = both.trial_correct.face{b}{t} +1; %correct specific to the trial
                    both.correct_keys.face{b}{t} = [both.correct_keys.face{b}{t} probe_resp(1)]; %note down frame correct key was pressed
                    
                    both.face_RT.secs{b}{t} = [both.face_RT.secs{b}{t} RT];
                    both.face_RT.frames{b}{t} = [both.face_RT.frames{b}{t} RT_id];
                    
                    both.face_resp_summary{b}{t}(i) = 1; %note it was correct
                else %if this was a house task
                    both.housecorrect = both.housecorrect + 1;
                    both.trial_correct.house{b}{t} = both.trial_correct.house{b}{t} +1; %correct specific to the trial
                    both.correct_keys.house{b}{t} = [both.correct_keys.house{b}{t} probe_resp(1)]; %note down frame correct key was pressed
                    
                    both.house_RT.secs{b}{t} = [both.house_RT.secs{b}{t} RT];
                    both.house_RT.frames{b}{t} = [both.house_RT.frames{b}{t} RT_id];
                    
                    both.house_resp_summary{b}{t}(i) = 1; %note it was correct
                end
                
                %remove that key from list so we don't count it twice
                trial_keys(resp_id(1)) = []; 
                trial_keySecs(resp_id(1)) = [];
                
            end
        end
        
        %if there are any leftover keys from the trial that didn't count towards a hit, then it's a false alarm
        both.trial_falsealarm{b}{t} = length(trial_keys);
        
        %update count for each stim
        both.face_count = both.face_count + both.changes{b}{t}/2;
        both.house_count = both.house_count + both.changes{b}{t}/2;
        
        both.FA_ratio{b}{t} = both.trial_falsealarm{b}{t}/length(both.keyi{b}{t}); % #false alarms/total responses
    end
end

both.faceaccuracy = both.facecorrect/both.face_count;
both.houseaccuracy = both.housecorrect/both.house_count;
both.mean_FAR = nanmean([both.FA_ratio{1}{:} both.FA_ratio{2}{:}]);
disp(['both, face accuracy: ', num2str(both.faceaccuracy)]);
disp(['both, house accuracy: ', num2str(both.houseaccuracy)]);
disp(['both, FAR: ', num2str(both.mean_FAR)]);
  
