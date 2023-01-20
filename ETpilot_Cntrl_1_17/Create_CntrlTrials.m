%% Creates valid trials for main experiment
% for control expt for RandDir_Trial
% currently coded to display trajectories for two stimuli within border
% one stim remains at fixation while other moves around with at least 75% overlap 

% will NOT include attend-both condition

% NOTE: this version replaces all 'house' instances with 'flower'

global Display Params

%% collect user inputted parameters

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
    q{n} = 'Practice? (0,1)';                     defaults{n} = '0';       n=n+1;
    q{n} = 'Eyetracking (0,1)';                   defaults{n} = '0';       n=n+1;
    q{n} = 'Room (0 = 425, 1 = mac, 2 = win)';    defaults{n} = '1';       

    answer = inputdlg(q,'Experimental Setup Information',1,defaults); 

    n=1;
    Params.info.SubIni = answer{n};                n=n+1;
    Params.info.SubNum = str2double(answer{n});    n=n+1;
    Params.info.FroBlo = str2double(answer{n});    n=n+1;
    Params.info.To_Blo = str2double(answer{n});    n=n+1;
    Params.info.Prac = str2double(answer{n});      n=n+1;
    Params.info.Eyetrk = str2double(answer{n});    n=n+1;
    Params.info.EyeRoom = str2double(answer{n});   



    if sum(isnan([Params.info.SubNum,Params.info.FroBlo,Params.info.To_Blo, Params.info.Prac, Params.info.Eyetrk, Params.info.EyeRoom])) ~=0
        errordlg('Please check parameters');
    else
        ExptStart = 1;
    end
end
    
Params.info.datetime = clock;
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
refresh = 1/60; 
Display.refresh = refresh;

ListenChar(2); %turn off keyboard input
HideCursor;

%trial parameters  
Params.randSeed_gen = rng('shuffle'); %can't run this on testing computer    
Params.trialLength = 10; % secs
Params.trialNum = 20;
Params.trialsGenerated = 100; 
Params.keyToResp = 'y';

% calculate pixels per degree 
Display.xCenter = Display.rect(3)/2; Display.yCenter = Display.rect(4)/2;
Display.v_dist = 60;   % viewing distance (cm)
Display.mon_width = 40.6;   % horizontal dimension of viewable screen (cm)
Display.ppd = pi * (Display.rect(3) - Display.rect(1)) / atan(Display.mon_width/Display.v_dist/2) / 360;

% dimensions of the box we want to restrict our movement
border_dist = 2; border = border_dist*Display.ppd;
border_rect = CenterRectOnPoint([0,0, Display.rect(3)-border,Display.rect(4)-border], Display.xCenter, Display.yCenter);

% stimulus: slight adjustment from main expt due to change in testing computer monitor
Params.stimSize = 11; %in degrees (300/ppd of old computer = 10.9526 rounded)

%% create image stimuli for main task

 
% specify which image stimuli to use here
s1 = 'face';
s2 = 'flower';

%starting face/flower difficulty (increased number for more difficulty) -
% doesn't matter what we set now, will change during practice
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

run OverlapImages_Cntrl
stimRect = stim_pos;
[stim_xCtr,stim_yCtr] = RectCenter(stimRect); % both stimuli will start at the center

s1_normal = eval(sprintf('gc_%stex', s1));
s2_normal = eval(sprintf('gc_%stex', s2));

%squeeze change
s1_chg = eval(sprintf('gc_%stex2', s1));
s2_chg = eval(sprintf('gc_%stex2', s2));

% Generate a huge set of possible trials and their trajectories 
% dictated by 'trialsGenerated' variable

% one stim will remain at fixation while other moves randomly
% same trajectory will be used twice

possTrials = []; % for moving stim
possTrials0 = []; % for stationary stim 
total_frames = round(Params.trialLength/Display.refresh);

% stationary stim, variables needed to generate trajectory for moving stim
stimRect_shifted = repmat(stimRect, [total_frames 1]);
xyMatrix = repmat([stim_xCtr;stim_yCtr], [1 total_frames+1]);

startGen = GetSecs();
for j = 1:Params.trialsGenerated
    
    % for moving stim
    run GenerateTrajectory3
    possTrials{j} = stimRect_shifted2;
   
    % for the stationary coordinates
    possTrials0{j} = stimRect_shifted;
    
    disp(['Trial done:', num2str(j)])
%     text = ['Trial done:', num2str(j)];
%     DrawFormattedText(Display.wPtr,text,'center','center',0);    
%     Screen('Flip',Display.wPtr);
end
endGen = GetSecs();
disp(['it took ', num2str(endGen-startGen), 'secs to generate ', num2str(Params.trialsGenerated), ' trajectories.'])

%make a separate 'nontrajectory' for when either stimulus doesn't need to move
stimRect_noshift = stimRect_shifted; 

% 
%% select low correlating trials
% 
corr_idx = calculate_stimCorr(possTrials0, possTrials);
pT = [];
for low = 1:length(corr_idx)
    pT{low} = possTrials{corr_idx(low)};
end

pT0 = possTrials0; % to have an analogous variable to pT and pT for stationary stimuli
%% filter out the ones that cross the border

%first stimulus
validTrialsIdx = CollectValidTraj(pT, border_rect);
available_Idx = validTrialsIdx;

%second stimulus
% validTrialsIdx2 = CollectValidTraj(pT2, border_rect);

% %make sure they overlap - EDIT: not needed because will always compare with stationary stim at fixation
% valTIdx = intersect(validTrialsIdx, validTrialsIdx2);

% shuffle order of trials - EDIT: will reuse the same trajectories for each stim
% trialsToRun = randsample(validTrialsIdx, length(validTrialsIdx));
% trialsToRun2 = randsample(validTrialsIdx2, length(validTrialsIdx2));
%% save relevant variables

filename = sprintf('%s/ETpilotTrials_%d_%d.mat', Params.dir, Params.sub.info, Params.info.datetime(2:3));
save(filename, 'border_rect', 'corr_idx', 'Display', 'face1', 'face_chg', 'flower1', 'flower_chg', 'overlaps', ...
    'Params', 'possTrials0', 'possTrials', 'pT0', 'pT', ...
    's1', 's2', 'stimRect', 'stimRect_noshift', 'stim_xCtr','stim_yCtr', ...
    'total_frames', 'validTrialsIdx', 'validTrialsIdx2', 'valTIdx', 'weights')


