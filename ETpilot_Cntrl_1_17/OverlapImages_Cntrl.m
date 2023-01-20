% NOTE: this version replaces all 'house' instances with 'flower'

global Display Params

%% image locations

img_dir = addpath([pwd '/images']);

%% load images

face1 = 'face1.png';
face_chg =  sprintf('face1_squeeze_0%d.png', img_values.face(face_diff));
flower1 = 'newflower.png';
flower_chg = sprintf('newflower_squeeze_0%d.png', img_values.flower(flower_diff));

%capture the alpha channel of the original picture - must be png image with
% transparent background
[~,~,face_alpha] = imread(face1); 
[~,~,face_chg_alpha] = imread(face_chg); 
[flower_norm,~,flower_alpha] = imread(flower1);
[flower_chg_norm,~,flower_chg_alpha] = imread(flower_chg); 


face_norm = rgb2gray(imread(face1));
face_chg_norm = rgb2gray(imread(face_chg));


% contrast normalize
% facech_norm = imadjust(face_gray, [], [0, .8]);
% face_norm = imadjust(face_gray, [], [.1, .9]);
% face_chg_norm = imadjust(face_chg_gray, [], [.1, .9]);
% 
% flowerch_norm = imadjust(flower_gray, [], [0, .8]);
% flower_norm = imadjust(flower_gray, [], [.1, .9]);
% flower_chg_norm = imadjust(flower_chg_gray, [], [.1, .9]);

% give normalized pictures transparent backgrounds
%for the regular pics
face_norm(:,:,2) = face_norm(:,:,1); 
face_norm(:,:,3) = face_norm(:,:,1); 
face_norm(:,:,4) = face_alpha; % add in the alpha channel

flower_norm(:,:,2) = flower_norm(:,:,1); 
flower_norm(:,:,3) = flower_norm(:,:,1); 
flower_norm(:,:,4) = flower_alpha; % add in the alpha channel

%and distorted
face_chg_norm(:,:,2) = face_chg_norm(:,:,1); 
face_chg_norm(:,:,3) = face_chg_norm(:,:,1); 
face_chg_norm(:,:,4) = face_chg_alpha; % add in the alpha channel

flower_chg_norm(:,:,2) = flower_chg_norm(:,:,1); 
flower_chg_norm(:,:,3) = flower_chg_norm(:,:,1); 
flower_chg_norm(:,:,4) = flower_chg_alpha; % add in the alpha channel
%% get sizes of images

% [face_iy, face_ix, face_iz]=size(face_gray);

dest_ix = Params.stimSize*Display.ppd; dest_iy = Params.stimSize*Display.ppd; %eventually convert to pixels*ppd
%% open screen

% Open a double buffered fullscreen window and draw a gray background
% and front and back buffers.

% Screen('Preference', 'SkipSyncTests', 1);
% screenNumber=min(Screen('Screens'));
% [wPtr, rect]=Screen('OpenWindow',screenNumber, 255/2, []);
% Priority(MaxPriority(wPtr));
% Screen(wPtr,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% xCenter = rect(3)/2; yCenter = rect(4)/2;
% 
% % Find the color values which correspond to white and black.
% black=BlackIndex(screenNumber);
% gray=GrayIndex(screenNumber); 
% 
%% make textures

% center image positions on screen
stim_pos = CenterRectOnPoint([0,0, dest_ix, dest_iy], Display.xCenter, Display.yCenter);

% opaque grayscale textures with transparent background
% g_facetex=Screen('MakeTexture', wPtr, face_gray);

gc_facetex=Screen('MakeTexture', Display.wPtr, face_norm);
gc_facetex2=Screen('MakeTexture', Display.wPtr, face_chg_norm);

gc_flowertex=Screen('MakeTexture', Display.wPtr, flower_norm);
gc_flowertex2=Screen('MakeTexture', Display.wPtr, flower_chg_norm);

%% draw texture with alpha overlay
% % 
% % 
% Screen('DrawTextures', wPtr, [gc_facetex, gc_flowertex],[], [stim_pos;stim_pos]', [], [], .5);
% Screen('Flip', wPtr);
% % Screen('FrameRect', wPtr, 0, fruit_pos, 1); %draw border 
% % Screen('FrameRect', wPtr, 0, flower_pos, 1); %draw border 
% WaitSecs(1);
% 
% Screen('DrawTextures', wPtr, [gc_facetex3, gc_flowertex3],[], [stim_pos;stim_pos]', [], [], .5);
% Screen('Flip', wPtr);
% WaitSecs(.2);
% 
% Screen('DrawTextures', wPtr, [gc_facetex, gc_flowertex],[], [stim_pos;stim_pos]', [], [], .5);
% Screen('Flip', wPtr);
% WaitSecs(3);
% sca;


