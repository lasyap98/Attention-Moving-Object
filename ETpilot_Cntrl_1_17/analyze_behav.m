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
                house_idxtask{b}{t}(end+1) = f;
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
%% calculate valid intervals for correct response
% by frame
for b = 1:length(house_idxstop)
    house.poss_resp_i{b} = {};
    for t = 1:length(house_idxstop{b})
        house.poss_resp_i{b}{t} = house_idxtask{b}{t};
        house.poss_resp_i{b}{t} = [];
        for i = 1:length(house_idxstop{b}{t})
%             house.poss_resp_i{b}{t} = [house.poss_resp_i{b}{t} house_idxstop{b}{t}(i)+1:house_idxstop{b}{t}(i) + resp_interval_frames(2)];
            house.poss_resp_i{b}{t} = [house.poss_resp_i{b}{t} house_idxstart{b}{t}(i)+resp_interval_frames(1):house_idxstop{b}{t}(i) + resp_interval_frames(2)];
            %can make resp starting from 100ms after stim starts, up til 1000ms after stim ends 
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


for b = 1:length(house.keyi)
    for t = 1:length(house.keyi{b})
        
        house.trial_correct{b}{t} = 0; %correct specific to the trial
        house.false_al{b}{t} = 0; %false alarms specific to the trial
        house.correct_keys{b}{t} = []; %which keys were correct specific to trial
        
        for i = 1:length(house.keyi{b}{t})
            
            if ismember(house.keyi{b}{t}(i), house.poss_resp_i{b}{t}) %correct response
                house.trial_correct{b}{t} = house.trial_correct{b}{t} +1; 
                house.correct_keys{b}{t} = [house.correct_keys{b}{t} house.keyi{b}{t}(i)];
            else %if not a member of valid responses, false alarm
                house.false_al{b}{t} = house.false_al{b}{t} +1;
            end
            
        end
        
        % changes/2 to get number of for only one stim
        if house.trial_correct{b}{t} > house.changes{b}{t}/2
            house.trial_correct{b}{t} = house.changes{b}{t}/2;
        end
        house.correct = house.correct + house.trial_correct{b}{t}; %add each trial's correct responses to total count
        
        if house.attend{b}{t}==2 %on trials where only the attended stimulus changes
            house.correct_sing = house.correct_sing + house.trial_correct{b}{t};
            house.sing_count = house.sing_count + house.changes{b}{t}/2;
        elseif house.attend{b}{t}==3 %on trials where both the stimuli changes
            house.correct_doub = house.correct_doub + house.trial_correct{b}{t};
            house.doub_count = house.doub_count + house.changes{b}{t}/2;
        end
        house.total_count = house.total_count + house.changes{b}{t}/2;
        house.FA_ratio{b}{t} = house.false_al{b}{t}/length(house.keyi{b}{t}); % #false alarms/total responses
    end
end

% for b = 1:length(house.keyi)
%     for t = 1:length(house.keyi{b})
%         
%         respTimes = house.keySecs_trial{b}{t}(house.keySecs_trial{b}{t}~=0);
%         house.correct_keys{b}{t} %the frame the correct key was pressed
%         assert(length(respTimes)==length(house.keyi{b}{t}), 'key_i not the same length as key_Secs!!')
%         
%         for i = 1:length(respTimes)
%             RT = respTimes-house.task.startTimes{b}{t}(i);
%             
%             if respTimes(i) >= house.possRT{b}{t}(i,1) && respTimes(i) <= house.possRT{b}{t}(i,2)
%                 house.correct = house.correct+1;
%                 
%                 if house.attend{b}{t}==2
%                     house.correct_sing = house.correct_sing+1; 
%                 elseif house.attend{b}{t}==3
%                     house.correct_doub = house.correct_doub+1; 
%                 end 
%                 
%             end
%         end
%         
%         if house.attend{b}{t}==2
%             house.sing_count = house.sing_count + house.changes{b}{t};
%         elseif house.attend{b}{t}==3
%             house.doub_count = house.doub_count + house.changes{b}{t};
%         end
%         house.total_count = house.total_count+ house.changes{b}{t};
%     end
% end
  

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

%% calculate valid intervals for correct response
for b = 1:length(face_idxstop)
    face.poss_resp_i{b} = {};
    for t = 1:length(face_idxstop{b})
%         face.poss_resp_i{b}{t} = face_idxtask{b}{t};
        face.poss_resp_i{b}{t} = [];
        for i = 1:length(face_idxstop{b}{t})
%             face.poss_resp_i{b}{t} = [face.poss_resp_i{b}{t} face_idxstop{b}{t}(i)+1:face_idxstop{b}{t}(i) + resp_interval_frames(2)];
            face.poss_resp_i{b}{t} =  [face.poss_resp_i{b}{t} face_idxstart{b}{t}(i)+resp_interval_frames(1) : face_idxstop{b}{t}(i) + resp_interval_frames(2)];
            %can make resp starting from 100ms after stim starts, up til 1000ms after stim ends
        end
        face.poss_resp_i{b}{t} = sort(face.poss_resp_i{b}{t});
    end
end
%% calculate accuracy - face
face.correct = 0;
face.correct_sing = 0;
face.correct_doub = 0;
face.doub_count = 0;
face.sing_count = 0;
face.total_count = 0;

for b = 1:length(face.keyi)
    for t = 1:length(face.keyi{b})
        
        face.trial_correct{b}{t} = 0;%correct specific to the trial
        face.false_al{b}{t} = 0; %false alarms specific to the trial
        face.correct_keys{b}{t} = []; %when correct keys were pressed specific to trial
        
        for i = 1:length(face.keyi{b}{t})
            
            if ismember(face.keyi{b}{t}(i), face.poss_resp_i{b}{t}) %correct response     
                face.trial_correct{b}{t} = face.trial_correct{b}{t} +1; 
                face.correct_keys{b}{t} = [face.correct_keys{b}{t} face.keyi{b}{t}(i)]; %note down frame correct key was pressed
            else %if not a member of valid responses, false alarm
                face.false_al{b}{t} = face.false_al{b}{t} +1; 
            end
            
        end
        
        % changes/2 to get number of for only one stim
        if face.trial_correct{b}{t} > face.changes{b}{t}/2
            face.trial_correct{b}{t} = face.changes{b}{t}/2;
        end
        face.correct = face.correct+face.trial_correct{b}{t}; %add each trial's correct responses to total count
        
        % separate accuracy based on whether trials included (1) only the attended stimulus changes or (2) where both the stimuli changes
        if face.attend{b}{t}==1
            face.correct_sing = face.correct_sing + face.trial_correct{b}{t}; 
            face.sing_count = face.sing_count + face.changes{b}{t}/2; 
        elseif face.attend{b}{t}==3
            face.correct_doub = face.correct_doub + face.trial_correct{b}{t}; 
            face.doub_count = face.doub_count + face.changes{b}{t}/2;
        end
        face.total_count = face.total_count+ face.changes{b}{t}/2;
        face.FA_ratio{b}{t} = face.false_al{b}{t}/length(face.keyi{b}{t}); % #false alarms/total responses
    end
end

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
        
%         evt_order{b}{t} = [];
%         for i = 1:both.changes{b}{t}
%             if both_faceidxstart{b}{t}(i) < both_houseidxstart{b}{t}(i)
%                 evt_order{b}{t} = [evt_order{b}{t} 1,2];
%             elseif both_faceidxstart{b}{t}(i) > both_houseidxstart{b}{t}(i)
%                 evt_order{b}{t} = [evt_order{b}{t} 2,1];
%             else
%                 evt_order{b}{t} = [evt_order{b}{t} randsample([1,2],2)];
%             end
%         end
    end
end


%% calculate valid intervals for correct response - both
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

%% accuracy - both
both.facecorrect = 0;
both.housecorrect = 0;
both.face_count = 0;
both.house_count = 0;

for b = 1:length(both_faceidxstart)
    both.FA_ratio_all{b} = [];
    for t = 1:length(both_faceidxstart{b})
        [all_start_idx, start_sort] = sort([both_faceidxstart{b}{t}, both_houseidxstart{b}{t}]);
        all_stop_idx = [both_faceidxstop{b}{t}, both_houseidxstop{b}{t}];
        all_stop_idx = all_stop_idx(start_sort);
        task_start = [ones([1,length(both_faceidxstart{b}{t})]) 2*ones([1,length(both_houseidxstart{b}{t})])];
        task_start = task_start(start_sort); %the order of events based on task
        
        trial_keys = both.keyi{b}{t};
        %initialize trial counts
        both.trial_correct.face{b}{t} = 0;
        both.trial_correct.house{b}{t} = 0;
        both.trial_falsealarm{b}{t} = 0;
        both.correct_keys.face{b}{t} = []; 
        both.correct_keys.house{b}{t} = []; 
        
        for i = 1:length(all_start_idx)
            % don't count the first .1 secs after stimulus starts 
            evt_resp_interval = all_start_idx(i)+resp_interval_frames(1):all_stop_idx(i)+resp_interval_frames(2);
            [probe_resp, resp_id] = intersect(trial_keys, evt_resp_interval);
            
            if ~isempty(resp_id) %if there is a valid response
                trial_keys(resp_id) = []; %remove that key from list so we don't count it twice
                if task_start(i) == 1 %if this was a face task
                    both.facecorrect = both.facecorrect + 1;
                    both.trial_correct.face{b}{t} = both.trial_correct.face{b}{t} +1; %correct specific to the trial
                    both.correct_keys.face{b}{t} = [both.correct_keys.face{b}{t} probe_resp]; %note down frame correct key was pressed
                else %if this was a house task
                    both.housecorrect = both.housecorrect + 1;
                    both.trial_correct.house{b}{t} = both.trial_correct.house{b}{t} +1; %correct specific to the trial
                    both.correct_keys.house{b}{t} = [both.correct_keys.house{b}{t} probe_resp]; %note down frame correct key was pressed
                end
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
  
