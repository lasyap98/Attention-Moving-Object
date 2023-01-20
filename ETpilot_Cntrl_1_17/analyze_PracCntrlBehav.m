function accuracy = analyze_PracCntrlBehav(sub, iBlock, attend, key_i, keyCode, taskframes, changes)

% 1/18: all instances of 'house' have been changed to flower

global Display Params pracOn lagOn

% RT should only count if they fall within a certain time-frame after change
% my average was around 500ms, set to 1 sec?
resp_interval = [.1,1]; %sec
resp_interval_frames = round(resp_interval/Display.refresh); %frames
change_frames = Params.change_frames;

% make sure we have the right block order
if pracOn
    block = sub.subjBlock(sub.subject,:);
elseif lagOn
    block = sub;
end

%% start stop idx of task

for t = 1:length(taskframes)
    idxtask{t} = [];
    idxstart{t} = [];
    idxstop{t} = {};
    for f = 1:length(taskframes{t})
        if taskframes{t}(block(iBlock),f)==1
            idxtask{t}(end+1) = f;
        end
    end
    len = length(idxtask{t});
    idxstart{t} = [idxtask{t}(1), idxtask{t}(change_frames+1), idxtask{t}((change_frames*2)+1)];
    idxstop{t} = [idxtask{t}(change_frames), idxtask{t}(change_frames*2), idxtask{t}(change_frames*3)];
    if (len/change_frames)==4
        idxstart{t}(end+1) = idxtask{t}((change_frames*3)+1);
        idxstop{t}(end+1) = idxtask{t}(change_frames*4);
    elseif (len/change_frames)==5
        idxstart{t}(end+1) = idxtask{t}((change_frames*3)+1);
        idxstart{t}(end+1) = idxtask{t}((change_frames*4)+1);
        
        idxstop{t}(end+1) = idxtask{t}(change_frames*4);
        idxstop{t}(end+1) = idxtask{t}(change_frames*5);
    end
    
    
end
%% calculate valid intervals for correct response

for t = 1:length(idxstop)
    poss_resp_i{t} = idxtask{t};
    for i = 1:length(idxstop{t})
        poss_resp_i{t} = [poss_resp_i{t} idxstop{t}(i)+1:idxstop{t}(i) + resp_interval_frames(2)];
    end
    poss_resp_i{t} = sort(poss_resp_i{t});
end

%% accuracy
evt_correct = 0;
correct_sing = 0;
correct_doub = 0;

doub_count = 0;
sing_count = 0;
total_count = 0;

for t = 1:length(key_i)
    trial_correct{t} = 0;%correct specific to the trial
    for i = 1:length(key_i{t})
        if ismember(key_i{t}(i), poss_resp_i{t})
            trial_correct{t} = trial_correct{t} +1; %correct specific to the trial
            if attend{t}==2
                correct_sing = correct_sing+1;
            elseif attend{t}==3
                correct_doub = correct_doub+1;
            end
        end
    end
    
    if trial_correct{t} > changes{t}/2
        trial_correct{t} = changes{t}/2;
    end
    %add number of correct responses in current trial to total count
    evt_correct = evt_correct+trial_correct{t};
    
    if attend{t}==2
        sing_count = sing_count + changes{t}/2;
    elseif attend{t}==3
        doub_count = doub_count + changes{t}/2;
    end
    total_count = total_count+ changes{t}/2;
end


accuracy = evt_correct/total_count;
correct_sing/sing_count;
correct_doub/doub_count;
% disp(['flower accuracy: ', num2str(flower.accuracy)]);

end
