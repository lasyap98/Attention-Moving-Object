if lagOn
%% sample trials to measure lag
% responses collected, only present one stim with behav task

global Display Params


%%

% create image stimuli
run OverlapImages_Cntrl
stimRect = stim_pos;
[stim_xCtr,stim_yCtr] = RectCenter(stimRect); % both stimuli will start at the center

% specify which image stimuli to use here
% s1 = 'face';
% s2 = 'flower';

s1_normal = eval(sprintf('gc_%stex', s1));
s2_normal = eval(sprintf('gc_%stex', s2));

%spherize change
s1_chg = eval(sprintf('gc_%stex2', s1));
s2_chg = eval(sprintf('gc_%stex2', s2));

%%

Lag.tperS = 5; %trials per stimulus

% shuffle order of trials
Lag_t2Run = randsample(validTrialsIdx, Lag.tperS);
Lag_t2Run2 = randsample(validTrialsIdx, Lag.tperS);
Lag.subjs = [1,2]; %face, flower

for b = 1:2
        
    for k = 1:Lag.tperS
        [TF{k}, CN{k}] = DetectChange(.2);
        Lag.attend{b}{1,k} = Lag.subjs(b); % 1 = attend face only, 2 = attend flower only
        
        if Lag.attend{b}{1,k} == 1 %face, single stim change
            Lag.attend{b}{2,k}  = 1; 
        elseif Lag.attend{b}{1,k} ==2 %flower, single stim change
            Lag.attend{b}{2,k} = 2;
        end
        
    end
    Lag.Task_block{b} = TF;
    Lag.Changes_block{b} = CN;
    
    Lag.t2Run{b} = Lag_t2Run;
    Lag.t2Run2{b} = Lag_t2Run2;
    
    clear TF CN
end

for s = 1:2
    if s==1
        stimulus = s1_normal;
        text = sprintf('New Block! \n\n Attend to the %s. \n\n Press ''%s'' each time you see the %s change. ', s1, Params.keyToResp, s1); %instructions
        s_name = s1;
    else
        stimulus = s2_normal;
        text = sprintf('New Block! \n\n Attend to the %s. Press ''%s'' each time you see the %s change.', s2, Params.keyToResp, s2); %instructions          
        s_name = s2;
    end
   
    
    DrawFormattedText(Display.wPtr,text,'center','center',0);    
    Screen('Flip',Display.wPtr);
    WaitSecs(.5);
    KbWait(); % wait for key press to begin trial

    ListenChar(2);
    HideCursor;
    
    
    if ET
        Lag.EyeData{s}.ELmsg_tm = [];
        Eyelink('message', sprintf('EVENT_LagStart_%s', s_name));
        Lag.EyeData{s}.ELmsg_tm = [Lag.EyeData{s}.ELmsg_tm GetSecs()];
    end

for iTrial = 1:Lag.tperS
    if ET
        Lag.EyeData{s}.mx{iTrial}=[];
        Lag.EyeData{s}.my{iTrial}=[];
        Lag.EyeData{s}.ma{iTrial}=[];
        Lag.EyeData{s}.time{iTrial}=[];
        Lag.EyeData{s}.FixDoneT{iTrial} = [];
    end
    
    Lag.currTrial_stim{iTrial} = pT{Lag_t2Run(iTrial)}; % retrieve the trajectory of a valid trial
    Lag.currTrial_stim2{iTrial} = pT2{Lag_t2Run2(iTrial)};
    
    if iTrial~= 1
        text = 'Press space to continue onto the next trial. Press ''q'' to quit.'; %instructions
        DrawFormattedText(Display.wPtr,text,'center','center',0);    
        Screen('Flip',Display.wPtr);
    end
    
    [~, keyCode, ~] = KbWait(); % wait for key press to begin trial
    if keyCode(KbName('q'))
        break
    end
    
    i = 1;
    
    %record responses
    Lag.keySecs_trial{s}{iTrial} = [];
    Lag.keyTime_trial{s}{iTrial} = [];
    Lag.keyCode_trial{s}{iTrial} = [];
    Lag.key_i{s}{iTrial} = [];

    KbQueueCreate();
    KbQueueStart();
    
    while i <= length(Lag.currTrial_stim{iTrial})

        Lag_EyeStart(iTrial) = GetSecs(); %time we start caring about eyetracking   
        currPosID = i;
        run FixCheck; %check eyetracker
        
        while  i==1
                stimPos = stimRect;
                stimPos2 = stimRect;
                %fixate for 500ms before proceeding
                Screen('DrawTexture', Display.wPtr, stimulus, [], [stimPos]', [], [], Params.alphaLevel);
%                 %fixate for 500ms before proceedingq
%                 Screen('FillOval', wPtr, [0 0 255] , stimPos); %draw blue circle, starting position
                Screen('FrameRect', Display.wPtr, 0, border_rect, Params.bordWidth); %draw borders
                Lag.trialStartT{s}(iTrial) = Screen('Flip',Display.wPtr);

                run FixCheck; %check eyetracker   
                if (GetSecs()-Lag_EyeStart(iTrial) > Params.fixDur) && (Fixated || ~ET)
                    i=i+1;   
                end
                Lag.flipTimes{s}{iTrial} = Lag.trialStartT{s}(iTrial);
        end
        
        if ET && i==2
            Eyelink('message', sprintf('EVENT_LagStartTrial%d_%s', iTrial, s_name));
            Lag.EyeData{s}.ELmsg_tm = [Lag.EyeData{s}.ELmsg_tm GetSecs()];
        end
        
        % start the trial trajectory
        
        switch Lag.attend{s}{2,iTrial}
            case 1 %change face only                    
                if Lag.Task_block{s}{iTrial}(1,i) == 1 %change task on for stim1
                    stimulus = s1_chg;
                else
                    stimulus = s1_normal;
                end
                stimPos = Lag.currTrial_stim{iTrial}(i,:);
            case 2 %change flower only
                if Lag.Task_block{s}{iTrial}(2,i) == 1 %change task on for stim2
                    stimulus = s2_chg;
                else
                    stimulus = s2_normal;
                end
                stimPos = Lag.currTrial_stim2{iTrial}(i,:);
        end
        
        % we only want to present 1 stimulus at a time   
        Screen('DrawTexture', Display.wPtr, stimulus, [], [stimPos]', [], [], Params.alphaLevel);
        Screen('FrameRect', Display.wPtr, 0, border_rect, Params.bordWidth); %draw border
        
        time = Screen('Flip', Display.wPtr);
        Lag.flipTimes{s}{iTrial} = [Lag.flipTimes{s}{iTrial} time];
        
        %record responses
        [pressed, firstPress] = KbQueueCheck();
        
        % KbName(firstPress)to generate a list of the keys for which the events occurred
        if pressed
            keyPress = KbName(find(firstPress)); 
            if ~iscell(keyPress) && strcmp(keyPress, Params.keyToResp) %|| strcmp(keyPress,'m')
                Lag.keyCode_trial{s}{iTrial} = [Lag.keyCode_trial{s}{iTrial} keyPress];
                Lag.key_i{s}{iTrial} = [Lag.key_i{s}{iTrial} i];
                Lag.keySecs_trial{s}{iTrial}  = [Lag.keySecs_trial{s}{iTrial} firstPress];
            end
            
        end
        
        i = i+1; %advance frame
    end
    
    ListenChar(0);
    KbQueueFlush();
    Lag.trialEndT{s}(iTrial) = Screen('Flip', Display.wPtr);
    Lag.trial_dur{s}(iTrial) = Lag.trialEndT{s}(iTrial)-Lag.trialStartT{s}(iTrial);
    Lag.IniFixDur{s}(iTrial) = Lag_InitialFixDur;
    
    if ET
        Eyelink('message', sprintf('EVENT_LagEndTrial%d_%s', iTrial, s_name));
        Lag.EyeData{s}.ELmsg_tm = [Lag.EyeData{s}.ELmsg_tm GetSecs()];
    end
    
end

if ET
    Eyelink('message', sprintf('EVENT_LagEnd_%s', s_name));
    Lag.EyeData{s}.ELmsg_tm = [Lag.EyeData{s}.ELmsg_tm GetSecs()];
end

Lag.accuracy{s} = analyze_PracCntrlBehav(Lag.subjs, s, Lag.attend{s}, Lag.key_i{s}, Lag.keyCode_trial{s}, ...
    Lag.Task_block{s}, Lag.Changes_block{s});

end
end
