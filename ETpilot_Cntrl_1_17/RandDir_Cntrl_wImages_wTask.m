%% Opens Screen, compiles valid trials, and displays on Screen
% control expt for RandDir_Trial
% currently coded to display trajectories for two stimuli within border
% one stim remains at fixation while other moves around with at least 75% overlap 

% NOTE: this version replaces all 'house' instances with 'flower'

global EyeData currPosID Display Params pracOn lagOn

%% collect user inputted parameters
%load the relevant trials
load('ET_pilotTrials_11_1.mat');

ExptStart = 0;

while ~ExptStart
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Dialog Box (can type in answers)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n=1;
    q{n} = 'Subject Initials';                    defaults{n} = 'test';    n=n+1;
    q{n} = 'Subject Number';                      defaults{n} = '999' ;    n=n+1;
    q{n} = 'From Block';                          defaults{n} = '1';       n=n+1;
    q{n} = 'To Block';                            defaults{n} = '4';       n=n+1;
    q{n} = 'Practice? (0,1)';                     defaults{n} = '1';       n=n+1;
    q{n} = 'Lag? (0,1)';                          defaults{n} = '1';       n=n+1;
    q{n} = 'Eyetracking? (0,1)';                  defaults{n} = '1';       n=n+1;
    q{n} = 'Room (0 = 425, 1 = mac, 2 = win)';    defaults{n} = '0';       

    answer = inputdlg(q,'Experimental Setup Information',1,defaults); 

    n=1;
    Params.info.SubIni = answer{n};                n=n+1;
    Params.info.SubNum = str2double(answer{n});    n=n+1;
    Params.info.FroBlo = str2double(answer{n});    n=n+1;
    Params.info.To_Blo = str2double(answer{n});    n=n+1;
    Params.info.Prac = str2double(answer{n});      n=n+1;
    Params.info.Lag = str2double(answer{n});       n=n+1;
    Params.info.Eyetrk = str2double(answer{n});    n=n+1;
    Params.info.EyeRoom = str2double(answer{n});   



    if sum(isnan([Params.info.SubNum,Params.info.FroBlo,Params.info.To_Blo, Params.info.Prac, Params.info.Lag, Params.info.Eyetrk, Params.info.EyeRoom])) ~=0
        errordlg('Please check parameters');
    else
        ExptStart = 1;
    end
end
    
Params.sub.info = sprintf('%s%d', Params.info.SubIni, Params.info.SubNum);
Params.dir = sprintf('data_cntrl/%s', Params.sub.info);
Params.info.datetime = clock;
mkdir(Params.dir);
%% set up room parameters

switch Params.info.EyeRoom
    case 0 %self
        Display.screen_coor = [];
        ETroom = 1;
        % eyetracking on (1) or off (0)
        ET =  Params.info.Eyetrk;
        downET = 1;
    case 1 %on mac running code for 425
        Display.screen_coor = [0 0 1024 768];
        ETroom = 0;
        ET = 0;
        downET = 0;
    case 2 %windows
        Display.screen_coor = [0 0 1920 1080];
        ETroom = 0;
        ET = 0;
        downET = 0;
end

%% open screen and set up trial parameters

Screen('Preference', 'SkipSyncTests', 1);
[Display.wPtr, Display.rect] = Screen('OpenWindow',min(Screen('Screens')), 255/2, Display.screen_coor);
Screen(Display.wPtr,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Priority(MaxPriority(Display.wPtr));
Screen('Flip', Display.wPtr);
refresh = Screen('GetFlipInterval',  Display.wPtr); % =.0167 = 60 Hz (for a 10s trial, sampled 600 times)
% refresh = 1/60;
Display.refresh = refresh;

ListenChar(2); %turn off keyboard input
HideCursor;

%trial parameters  
Params.randSeed = rng('shuffle'); 
Params.trialLength = 10; % secs
Params.trialNum = 20;
Params.trialsGenerated = 100; 
Params.keyToResp = 'y';
Params.fixDur = .5;
Params.bordWidth = 5;
Params.alphaLevel = .5;

% calculate pixels per degree 
Display.xCenter = Display.rect(3)/2; Display.yCenter = Display.rect(4)/2;
Display.v_dist = 60;   % viewing distance (cm)
Display.mon_width = 40.6;   % horizontal dimension of viewable screen (cm)
Display.ppd = pi * (Display.rect(3) - Display.rect(1)) / atan(Display.mon_width/Display.v_dist/2) / 360;


%%

if ET
    EyelinkSetup(0);
end
Fixated = 0;

 
% specify which image stimuli to use here
s1 = 'face';
s2 = 'flower';
%% run single stim trials for measure ET lag 

%starting face/flower difficulty (increased number for more difficulty)
face_diff = 1; %1-8
flower_diff = 1; %1-8

if Params.info.Prac
    faceOK = 0; %have NOT finalized the difficulty
    flowerOK = 0;
else 
    faceOK = 1; %assuming practice has already been completed (or is not needed)
    flowerOK = 1;
end

if (~faceOK || ~flowerOK)
    %starting level fo naive
    img_values.face = [15, 10, 6, 3]; 
    img_values.flower =  [15, 10, 6, 3]; 
else
    img_values.face = [15, 12, 10, 8, 7, 6, 5, 4, 3]; 
    img_values.flower = [15, 12, 10, 8, 7, 6, 5, 4, 3];
end

pracOn = 0; %don't start Practice
lagOn = Params.info.Lag; % start Lag
run Cntrl_measureLag
lagOn = 0; %turn off Lag

%% practice

if Params.info.Prac
    pracOn = 1;
    run PracticeCntrlTrials
    pracOn = 0;
    sca;
    %this will quit the screen
end

%% plot practice performance

% calculate full range of performance
Titrate.full_fitValues.face = polyval(Titrate.fitParams.face, img_values.face);
Titrate.full_fitValues.flower = polyval(Titrate.fitParams.flower, img_values.flower);

%plot face
plot(Titrate.face_diffs_all, [Titrate.faceAcc{:}], 'o')
hold on
plot(Titrate.face_diffs_all, Titrate.fitValues.face)
hold on
plot(img_values.face, Titrate.full_fitValues.face)
hold on
%plot flower
plot(Titrate.flower_diffs_all, [Titrate.flowerAcc{:}], 'o')
hold on
plot(Titrate.flower_diffs_all, Titrate.fitValues.flower)
hold on
plot(img_values.flower, Titrate.full_fitValues.flower)

desired_perctge = .75; % select the level closest to 75% accuracy - maybe don't consider if passes 75%?
[~, flower_diff_i] = min(abs(Titrate.full_fitValues.flower-desired_perctge)); %find the index closest to difficulty

[~, face_diff_i] = min(abs(Titrate.full_fitValues.face-desired_perctge)); %find the index closest to difficulty

% Use the same level for both stim or different levels
face_diff = img_values.face(face_diff_i);
flower_diff = img_values.flower(flower_diff_i);


%% create image stimuli for main task

%degrees, same for x and y
run OverlapImages_Cntrl
stimRect = stim_pos;
[stim_xCtr,stim_yCtr] = RectCenter(stimRect); % both stimuli will start at the center

s1_normal = eval(sprintf('gc_%stex', s1));
s2_normal = eval(sprintf('gc_%stex', s2));

%squeeze change
s1_chg = eval(sprintf('gc_%stex2', s1));
s2_chg = eval(sprintf('gc_%stex2', s2));

%% create main task
Params.sub.subjBlock(1,:) = [1, 2, 1, 2]; %face, flower, both
Params.sub.subjBlock(2,:) = [2, 1, 2, 1]; %flower, both, face

Params.sub.subject = mod(Params.info.SubNum, 2)+1;

for b = Params.info.FroBlo:Params.info.To_Blo
      
    % half of blocks should change just 1 stim and half should change both
    flip_chg_count = randsample([ones(1,Params.trialNum/2) zeros(1,Params.trialNum/2)], Params.trialNum); % 0 = both stim change, 1 = single stim
    
    for t = 1:Params.trialNum
        [TF{t}, CN{t}] = DetectChange(.2);
        
        attend{b}{1,t} = Params.sub.subjBlock(Params.sub.subject,b); % 1 = attend face only, 2 = attend flower only
        
        if flip_chg_count(t) && attend{b}{1,t} ==1 %face, single stim change
            attend{b}{2,t}  = 1; 
        elseif  flip_chg_count(t) && attend{b}{1,t} ==2 %flower, single stim change
            attend{b}{2,t} = 2;
        else %face, both; flower, both; both, both
            attend{b}{2,t} = 3; % 1 = change face only, 2 = change flower only, 3 = change both
        end
        
    end
    Task_block{b} = TF;
    Changes_block{b}= CN; 
end

%% start block
ListenChar(2); %turn off keyboard input


for iBlock = Params.info.FroBlo:Params.info.To_Blo
switch attend{iBlock}{1,1} %specify what image should be attended
    case 1
        text = sprintf('New Block! \n\n Attend to the %s. \n\n Press ''%s'' each time you see the %s change. \n\n Press space to continue.', s1, Params.keyToResp,  s1); %instructions
    case 2 
        text = sprintf('New Block! \n\n Attend to the %s. Press ''%s'' each time you see the %s change. \n\n Press space to continue.',  s2, Params.keyToResp, s2); %instructions          
end

DrawFormattedText(Display.wPtr,text,'center','center',0);    
Screen('Flip',Display.wPtr);
WaitSecs(.5);

% KbWait(); % wait for key press to begin trial
[~, keyCode] = KbWait(-1);
while keyCode(KbName('space'))==0
    [~, keyCode] = KbWait(-1);
end 

%start the block - may have to change this to a while loop?
Block.startT{iBlock} = GetSecs();
task_frames = Task_block{iBlock};

% shuffle order of trials at the start of each odd block: each stim needs
% to do the same trajectory before reshuffling
if mod(iBlock,2) == 1 % complementary attend condition (face after flower, or flower after face)
%     trialsToRun = randsample(validTrialsIdx, length(validTrialsIdx)); %shuffle the list
    t2Run = randsample(available_Idx, Params.trialNum); %select random sample
end

for iTrial = 1:Params.trialNum
    if ET
        EyeData.mx{iTrial}=[];
        EyeData.my{iTrial}=[];
        EyeData.ma{iTrial}=[];
        EyeData.time{iTrial}=[];
        EyeData.FixDoneT{iTrial} = [];
        EyeData.ELmsg_tm{iTrial} = [];
    end
    
    switch Block.attend{iBlock}{1,1} %specify what image(s) the eyes should be tracking
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
    end
    
    WaitSecs(.5);
    [~, keyCode, ~] = KbWait(); % wait for key press to begin trial
    if keyCode(KbName('q'))
        break
    end
    
    % initialize frame count
    i = 1;
    
    % to record responses
    keySecs_trial{iTrial} = [];
    keyTime_trial{iTrial} = [];
    keyCode_trial{iTrial} = [];
    key_i{iTrial} = [];
    
    %start collecting response
    KbQueueCreate();
    KbQueueStart();
    
    while i <= length(currTrial_stim{iTrial})

        if ET
            Eyelink('message', sprintf('EVENT_Block%dTrial%d_Start', iBlock, iTrial));
            EyeData.ELmsg_tm{iTrial} = [EyeData.ELmsg_tm{iTrial} GetSecs()]; %timestamp when we send Eyelink a message
        end
        
        EyeStart(iTrial) = GetSecs(); %time we start caring about eyetracking   
        currPosID = i;
        run FixCheck; %check eyetracker
        
        while  i==1
                stimPos = stimRect;
                stimPos2 = stimRect;
                
                Screen('DrawTextures', Display.wPtr, [s2_normal, s1_normal], [], [stimPos2;stimPos]', [], [], Params.alphaLevel);
                Screen('FrameRect', Display.wPtr, 0, border_rect, Params.bordWidth); %draw borders
                trialStartT(iTrial) = Screen('Flip',Display.wPtr);
                
                if ET
                    Eyelink('message', sprintf('EVENT_Block%dTrial%d_Fix', iBlock, iTrial));
                    EyeData.ELmsg_tm{iTrial} = [EyeData.ELmsg_tm{iTrial} GetSecs()]; %timestamp when we send Eyelink a message
                end

                run FixCheck; %check eyetracker   
                if (GetSecs()-EyeStart(iTrial) > Params.fixDur) && (Fixated || ~ET) %fixate for 500ms before proceeding
                    i=i+1;   
                end
                flipTimes{iTrial} = trialStartT(iTrial);
        end
        
        % start the trial trajectory
        stimPos = currTrial_stim{iTrial}(i,:);
        stimPos2 = currTrial_stim2{iTrial}(i,:);
        
        switch attend{iBlock}{2,iTrial}
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
        
        %record task event in ET file
        if ET 
            if attend{iBlock}{1,iTrial}==1 && stimulus1 == s1_chg %only let ET know if attend-face and face task on
                Eyelink('message', 'EVENT_TaskFaceOn');
                EyeData.ELmsg_tm{iTrial} = [EyeData.ELmsg_tm{iTrial} GetSecs()]; %timestamp when we send Eyelink a message
            elseif attend{iBlock}{1,iTrial}==2 && stimulus2 == s2_chg %only let ET know if attend-flower and flower task on
                Eyelink('message', 'EVENT_TaskHouseOn');
                EyeData.ELmsg_tm{iTrial} = [EyeData.ELmsg_tm{iTrial} GetSecs()]; %timestamp when we send Eyelink a message
            elseif (attend{iBlock}{1,iTrial}==0 && stimulus2 == s2_chg) %let ET know if attend-both and flower task on
                Eyelink('message', 'EVENT_TaskHouseOn');
                EyeData.ELmsg_tm{iTrial} = [EyeData.ELmsg_tm{iTrial} GetSecs()]; %timestamp when we send Eyelink a message
            elseif(attend{iBlock}{1,iTrial}==0 && stimulus1 == s1_chg) %let ET know if attend-both and face task on
                Eyelink('message', 'EVENT_TaskFaceOn');
                EyeData.ELmsg_tm{iTrial} = [EyeData.ELmsg_tm{iTrial} GetSecs()]; %timestamp when we send Eyelink a message
            end
        end
        
        %record responses
        [pressed, firstPress] = KbQueueCheck();
        
        % KbName(firstPress)to generate a list of the keys for which the events occurred
        if pressed
            keyPress = KbName(find(firstPress)); 
            if ~iscell(keyPress) && strcmp(keyPress,Params.keyToResp) %|| strcmp(keyPress,'m')
                %figure out what to do when participant presses space, bc it breaks code
                keyCode_trial{iTrial} = [keyCode_trial{iTrial} keyPress];
                key_i{iTrial} = [key_i{iTrial} i];
                keySecs_trial{iTrial}  = [keySecs_trial{iTrial} firstPress]; 
            end
            
        end

        i = i+1;
    end
        ListenChar(0); 
        KbQueueFlush();
        trialEndT(iTrial) = Screen('Flip', Display.wPtr);
        Block.trialDur{iBlock}(iTrial) = trialEndT(iTrial)-trialStartT(iTrial);
        
        if ET
            Eyelink('message', sprintf('EVENT_Block%dTrial%d_End', iBlock, iTrial));
            EyeData.ELmsg_tm{iTrial} = [EyeData.ELmsg_tm{iTrial} GetSecs()]; %timestamp when we send Eyelink a message
        end
        
end
% disp(['Trials lasted an average of ', num2str(mean(trial_durations(trial_durations~=0))), 'secs'])

%Block variables
%timings
Block.endT{iBlock} = GetSecs();
Block.IniFixDur{iBlock} = InitialFixDur;
Block.trialStartT{iBlock} = trialStartT;
Block.trialEndT{iBlock} = trialEndT;
Block.flipTimes{iBlock} = flipTimes;
%stimuli
Block.t2Run{iBlock} = t2Run;
Block.currTstim{iBlock} = currTrial_stim;
Block.currTstim2{iBlock} = currTrial_stim2;
%responses/ET
Block.keyCode_trial{iBlock}  = keyCode_trial;
Block.key_i{iBlock}  = key_i;
Block.keySecs_trial{iBlock}  = keySecs_trial;
Block.EyeData{iBlock} = EyeData;

available_Idx = setdiff(validTrialsIdx, t2Run); %remove the trial idx we just presented so to prevent repeat in next block
end


sca;
Screen('CloseAll'); % close all open onscreen and offscreen windows and textures, movies and videosources - release nearly all resources

ShowCursor();

if ET
    Eyelink('StopRecording');   
    Eyelink('CloseFile');
    
    %download edf file
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

% run analyze_behav %EITHER UPDATE CODE OR REMOVE

if ETroom %we only really care about the data if we run it in the ET room, not on my laptop
    filename = sprintf('%s/ETpilot_%s_%d_%d_%d_%d.mat', Params.dir, Params.sub.info, Params.info.datetime(2:5));
    save(filename)
end
