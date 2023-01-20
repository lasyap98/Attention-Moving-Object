%% Opens Screen, compiles valid trials, and displays on Screen
% currently coded to display trajectories for two stimuli within border
% with at least 75% overlap throughout trajectories  

global EyeData currPosID Display Params
%% begin practice

if pracOn
        
    text = 'Starting Practice Trials...';
    DrawFormattedText(Display.wPtr,text,'center','center',0);    
    Screen('Flip', Display.wPtr);

%     if ET
%         EyelinkSetup(0); 
%     end

    %trial parameters     
    Practice.trialNum = 5;
    Practice.blockNum = 2; %1 block per condition
    
    % create image stimuli
    run OverlapImages_Ctrl
    stimRect = stim_pos;
    [stim_xCtr,stim_yCtr] = RectCenter(stimRect); % both stimuli will start at the center

    s1_normal = eval(sprintf('gc_%stex', s1));
    s2_normal = eval(sprintf('gc_%stex', s2));

    %spherize change
    s1_chg = eval(sprintf('gc_%stex2', s1));
    s2_chg = eval(sprintf('gc_%stex2', s2));
    
    %% create task
    sub.subjBlock(1,:) = [1, 2]; %face, flower
    sub.subjBlock(2,:) = [2, 1]; %flower, face

    sub.subject = mod(Params.info.SubNum, 2)+1;

for b = 1:Practice.blockNum 
      
    % half of blocks should change just 1 stim and half should change both
    flip_chg_count = randsample([ones(1,3) zeros(1,2)], Practice.trialNum); % 0 = both stim change, 1 = single stim
    
    for k = 1:Practice.trialNum
        [TF{k}, CN{k}] = DetectChange(.2);
         
        Practice.attend{b}{1,k} = sub.subjBlock(sub.subject,b); % 1 = attend face only, 2 = attend flower only
        
        if flip_chg_count(k) && Practice.attend{b}{1,k} ==1 %face, single stim change
            Practice.attend{b}{2,k}  = 1; 
        elseif  flip_chg_count(k) && Practice.attend{b}{1,k} ==2 %flower, single stim change
            Practice.attend{b}{2,k} = 2;
        else %face, both; flower, both
            Practice.attend{b}{2,k} = 3; % 1 = change face only, 2 = change flower only, 3 = change both
        end
        
    end
    
    Practice.Task_block{b} = TF;
    Practice.Changes_block{b}= CN; 
    
    clear TF CN
    
    
end
%%

Prac_run = 1;
EyeData{Prac_run}.ELmsg_tm = [];
while (~faceOK || ~flowerOK)
    %% start block (Practice run)

    if ET
        Eyelink('message', sprintf('EVENT_PracRun%d_Start', Prac_run));
        EyeData{Prac_run}.ELmsg_tm = [EyeData{Prac_run}.ELmsg_tm GetSecs()];
    end

    % randomly select trials to run from possible trajectories
%     trialsToRun = randsample(validTrialsIdx, length(validTrialsIdx)); %shuffle the list
    t2Run = randsample(validTrialsIdx, Practice.trialNum); %select random sample
    
for iBlock = 1:Practice.blockNum
    switch Practice.attend{iBlock}{1,1} %specify what image(s) the eyes should be tracking
        case 1
            text = sprintf('New Practice Block! \n\n Attend to the %s. \n\n Press ''%s'' each time you see the %s change. \n\n Press space to continue.', s1, Params.keyToResp, s1); %instructions
        case 2
            text = sprintf('New Practice Block! \n\n Attend to the %s. Press ''%s'' each time you see the %s change. \n\n Press space to continue.', s2, Params.keyToResp, s2); %instructions
    end

    DrawFormattedText(Display.wPtr,text,'center','center',0);    
    Screen('Flip',Display.wPtr);
    WaitSecs(.5);
    % KbWait(); % wait for key press to begin trial
    
    [~, keyCode] = KbWait(-1);
    while keyCode(KbName('space'))==0
        [~, keyCode] = KbWait(-1);
    end 

    %start the trials 
    trial_durations = zeros([1,Practice.trialNum]);
    task_frames = Practice.Task_block{iBlock};
    
    %shuffle the selected list of trials to run before beginning block
%     trialsToRun = randsample(trialsToRun, length(trialsToRun)); 
    
    for iTrial = 1:Practice.trialNum
        if ET
            EyeData.mx{iTrial}=[];
            EyeData.my{iTrial}=[];
            EyeData.ma{iTrial}=[];
            EyeData.time{iTrial}=[];
            EyeData.FixDoneT{iTrial} = [];
        end
        
        switch Practice.attend{iBlock}{1,1} %specify what image(s) the eyes should be tracking
            case 1 %attend face
                currTrial_stim{iTrial} = pT{t2Run(iTrial)}; % face gets moving trajectory
                currTrial_stim2{iTrial} = pT0{t2Run(iTrial)}; % flower gets stationary traj
            case 2 %attend flower
                currTrial_stim{iTrial} = pT0{t2Run(iTrial)}; % face gets stationary traj
                currTrial_stim2{iTrial} = pT{t2Run(iTrial)}; % flower gets moving trajectory
        end
        
        if iTrial~= 1
            text = 'Press space to continue onto the next trial. Press ''q'' to quit.'; %instructions
            DrawFormattedText(Display.wPtr,text,'center','center',0);    
            Screen('Flip',Display.wPtr);
        
 
            WaitSecs(.5);
            [~, keyCode, ~] = KbWait(); % wait for key press to begin trial
            if keyCode(KbName('q'))
                break
            end
        end
        i = 1;

        %record responses
        keySecs_trial{iTrial} = [];
        keyCode_trial{iTrial} = [];
        key_i{iTrial} = [];

        KbQueueCreate();
        KbQueueStart();

        while i <= length(currTrial_stim{iTrial})
            
            EyeStart(iTrial) = GetSecs(); %time we start caring about eyetracking   
            currPosID = i;
            run FixCheck; %check eyetracker

            while  i==1
                    stimPos = stimRect;
                    stimPos2 = stimRect;
                    %fixate for 500ms before proceeding
                    Screen('DrawTextures', Display.wPtr, [s2_normal, s1_normal], [], [stimPos2;stimPos]', [], [], Params.alphaLevel);
    %                 %fixate for 500ms before proceedingq
    %                 Screen('FillOval', Display.wPtr, [0 0 255] , stimPos); %draw blue circle, starting position
                    Screen('FrameRect', Display.wPtr, 0, border_rect, Params.bordWidth); %draw borders
                    trialStartT = Screen('Flip',Display.wPtr);

                    run FixCheck; %check eyetracker   
                    if (GetSecs()-EyeStart(iTrial) > Params.fixDur) && (Fixated || ~ET)
                        i=i+1;   
                    end
                    flipTimes{iTrial} = trialStartT;
            end
            
            if ET && i==2
                Eyelink('message', sprintf('EVENT_PracRun%d_StartTrial%d_%s', Prac_run, iTrial, s1));
                EyeData{Prac_run}.ELmsg_tm = [EyeData{Prac_run}.ELmsg_tm GetSecs()];
            end

            % start the trial trajectory
            stimPos = currTrial_stim{iTrial}(i,:); %face
            stimPos2 = currTrial_stim2{iTrial}(i,:); %flower

            switch Practice.attend{iBlock}{2,iTrial}
                case 3 %change both
                    if task_frames{iTrial}(1,i) == 1 %change task on for stim1
                        stimulus1 = s1_chg;
                    else
                        stimulus1 = s1_normal;
                    end
                    if task_frames{iTrial}(2,i) == 1 %change task on for stim2
                        stimulus2 = s2_chg;
                    else
                        stimulus2 = s2_normal;
                    end
                case 1 %change face only                    
                    if task_frames{iTrial}(1,i) == 1 %change task on for stim1
                        stimulus1 = s1_chg;
                    else
                        stimulus1 = s1_normal;
                    end
                    stimulus2 = s2_normal;
                case 2 %change flower only
                    if task_frames{iTrial}(2,i) == 1 %change task on for stim2
                        stimulus2 = s2_chg;
                    else
                        stimulus2 = s2_normal;
                    end
                    stimulus1 = s1_normal;
            end

            Screen('DrawTextures', Display.wPtr, [stimulus2, stimulus1], [], [stimPos2;stimPos]', [], [], Params.alphaLevel);
            Screen('FrameRect', Display.wPtr, 0, border_rect, Params.bordWidth); %draw border

            time = Screen('Flip', Display.wPtr);
            flipTimes{iTrial} = [flipTimes{iTrial} time];

            %record responses
            [pressed, firstPress] = KbQueueCheck();

            % KbName(firstPress)to generate a list of the keys for which the events occurred
            if pressed
                keyPress = KbName(find(firstPress)); 
                if ~iscell(keyPress) && strcmp(keyPress,Params.keyToResp) 
                    keyCode_trial{iTrial} = [keyCode_trial{iTrial} keyPress];
                    key_i{iTrial} = [key_i{iTrial} i];
                    keySecs_trial{iTrial}  = [keySecs_trial{iTrial} firstPress];
                end

            end

            i = i+1; %advance to next frame
        end
            ListenChar(0); 
            KbQueueFlush();
            trialEndT = Screen('Flip', Display.wPtr);
            trial_durations(iTrial) = trialEndT()-trialStartT;
            
            if ET
                Eyelink('message', sprintf('EVENT_PracRun%d_EndTrial%d_%s', Prac_run, iTrial, s1));
                EyeData{Prac_run}.ELmsg_tm = [EyeData{Prac_run}.ELmsg_tm GetSecs()];
            end

    end
    disp(['Trials lasted an average of ', num2str(mean(trial_durations(trial_durations~=0))), 'secs'])

    if ET
        Eyelink('message', sprintf('EVENT_PracRun%d_End', Prac_run));
        EyeData{Prac_run}.ELmsg_tm = [EyeData{Prac_run}.ELmsg_tm GetSecs()];
    end
    
    %Block variables
    Practice.EyeData{iBlock} = EyeData;

    Practice.keyCode_trial{iBlock}  = keyCode_trial;
    Practice.key_i{iBlock}  = key_i;
    Practice.keySecs_trial{iBlock}  = keySecs_trial;

    Practice.flipTimes{iBlock} = flipTimes;
    Practice.trialDur{iBlock} = trial_durations;
    Practice.t2Run{iBlock} = trialsToRun;
    
    Practice.currTstim{iBlock} = currTrial_stim;
    Practice.currTstim2{iBlock} = currTrial_stim2;
    
    % calculate accuracy
    Practice.accuracy{iBlock} = analyze_PracCntrlBehav(sub, iBlock, Practice.attend{iBlock}, Practice.key_i{iBlock}, Practice.keyCode_trial{iBlock},...
        Practice.Task_block{iBlock}, Practice.Changes_block{iBlock});
    
    % assess accuracy and determine whether to change difficulty or not
    if Practice.attend{iBlock}{1,k} == 1
        Practice.face_acc = Practice.accuracy{iBlock};
        text = sprintf('Accuracy: %2.f%%', Practice.accuracy{iBlock}*100);
    elseif Practice.attend{iBlock}{1,k} == 2
        Practice.flower_acc = Practice.accuracy{iBlock};
        text = sprintf('Accuracy: %2.f%%', Practice.accuracy{iBlock}*100);
    end
    
    DrawFormattedText(Display.wPtr,text,'center','center',0);    
    Screen('Flip', Display.wPtr);
    WaitSecs(2);
    
end %iBlock

    % note down values for titration
    Titrate.faceAcc{Prac_run} = Practice.face_acc;
    Titrate.flowerAcc{Prac_run} = Practice.flower_acc;
    Titrate.face_diffs_all = img_values.face;
    Titrate.flower_diffs_all = img_values.flower;
    
    
    if  (face_diff ~= length(img_values.face)) % && (Practice.face_acc >= .8)
        face_diff =  face_diff + 1; %increase the difficulty
    elseif (face_diff == length(img_values.face)) % && (Practice.face_acc >= .8)
        faceOK = 1;
    end
       
    if  (flower_diff ~= length(img_values.flower)) % && (Practice.flower_acc >= .8) 
        flower_diff =  flower_diff + 1; %increase the difficulty
    elseif (flower_diff == length(img_values.flower)) % && (Practice.flower_acc >= .8)
        flowerOK = 1;
    end
    
    if ETroom %we only really care about the data if we run it in the ET room, not on my laptop
        filename = sprintf('%s/ETpilot_PracticeRun%d_%s_%d_%d_%d_%d.mat', Params.dir, Prac_run, Params.sub.info, Params.info.datetime(2:5));
        save(filename, 'Params', 'Practice', 'Lag', 'face_chg', 'flower_chg', 'Titrate')
    end
    
    %update images
    run OverlapImages_Ctrl
    stimRect = stim_pos;
    [stim_xCtr,stim_yCtr] = RectCenter(stimRect); % both stimuli will start at the center

    s1_normal = eval(sprintf('gc_%stex', s1));
    s2_normal = eval(sprintf('gc_%stex', s2));

    %spherize change
    s1_chg = eval(sprintf('gc_%stex2', s1));
    s2_chg = eval(sprintf('gc_%stex2', s2));
    
    %run the next level of practice
    Prac_run = Prac_run+1;
    
    
    
end

%% fit polynomial to practice performance
%fit polynomial to (number of runs - 1) terms

Titrate.fitParams.flower = polyfit(Titrate.flower_diffs_all, [Titrate.flowerAcc{:}], Prac_run-1); 
Titrate.fitValues.flower = polyval(Titrate.fitParams.flower, Titrate.flower_diffs_all);

Titrate.fitParams.face = polyfit(Titrate.face_diffs_all, [Titrate.faceAcc{:}], Prac_run-1);
Titrate.fitValues.face = polyval(Titrate.fitParams.face, Titrate.face_diffs_all);


%%
ShowCursor();

if ETroom %we only really care about the data if we run it in the ET room, not on my laptop    
    filename = sprintf('%s/ETpilot_EndPractice_%s_%d_%d_%d_%d.mat', Params.dir, Params.sub.info, Params.info.datetime(2:5));
    save(filename)
end

Screen('CloseAll');

if ET
    Eyelink('StopRecording');   
    Eyelink('CloseFile');
    
    %download edf file - this one will be separate to the ET data for the main task
    if downET
        Eyelink('Command', 'set_idle_mode');
    try
        fprintf('Receiving data file ''%s''\n',  EyeData.edfFile);
        status = Eyelink('ReceiveFile');
        
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(EyeData.edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n',  EyeData.edfFile, Params.dir );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n',  EyeData.edfFile);
    end
    end
    Eyelink('Shutdown');
end

end %while pracOn
