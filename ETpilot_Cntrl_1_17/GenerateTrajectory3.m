global Params Display
%% Generates random trajectory for ONE trial based on traj of first stimulus
% to ensure atleast 75% overlap with first stimulus (correct usage of bboxOverlap)

%% weights calculated with gaussian dist
x = linspace(-3,3,64);
sd = 1;
y = exp(-0.5*(x/sd).^2);
weights = (exp(-0.5*(x/sd).^2))/4; % to be assigned to 'n=length' previous directions
%% mean theta will need to be weighted average of n previous directions
% collect past n directions to include in weighted (gaussian) average 

speed = 1; % in pixels: the bigger the stimulus, the faster it will need to go   
total_frames = round(Params.trialLength/Display.refresh); % will be based on refresh rate
randn2 = rand([1,total_frames+length(weights)-1]); %how many random thetas do we want
randTheta2 = 2*pi*randn2; %create the random thetas


%% convert to delta-x and delta-y values

s2_delta_x2 = speed*cos(randTheta2);
s2_delta_y2 = speed*sin(randTheta2);

%% create matrix of weights, should have full table adjusted for each frame
% the greatest weight should be assigned to the most recent frame and shift accordingly

weightMatrix2 = zeros(length(randTheta2));

for i = 1:length(randTheta2)
    for k = 0:length(weights)-1
        if (i+k) > length(randTheta2)
            break
        end
        weightMatrix2(i+k,i) = weights(k+1); 
    end
end        
    
    
%% multiply matrices to get weighted delta-x and delta-y values
% 
weighted_x2 = s2_delta_x2*weightMatrix2'; %transposed matrix to put oldest on top and help calculate in sequential order
weighted_y2 = s2_delta_y2*weightMatrix2';

% only keep the ones with full weights calculated
weighted_x_trimmed2 = weighted_x2((length(weights)):end);
weighted_y_trimmed2 = weighted_y2((length(weights)):end);

%% update position matrix

xyMatrix2 = zeros([2,length(weighted_x_trimmed2)]); % (1,:) = x, (2,:) = y
xyMatrix2(:,1) = [stim_xCtr;stim_yCtr]; % starting position

stimRect_shifted2 = OffsetRect(stimRect, weighted_x2(1), weighted_y2(1)); %stimRect = center position of screen
[xyMatrix2(1,2),xyMatrix2(2,2)] = RectCenter(stimRect_shifted(1,:)); % xyMatrix is 1 ahead of stimRect_shifted bc it includes starting pos

%need to rearrange rect coordinates to [x(ctr), y(ctr), width, height] = 'bounding box'
box2 = [xyMatrix2(1,2), xyMatrix2(2,2), ...
    RectWidth(stimRect_shifted2(1,:)), RectHeight(stimRect_shifted2(1,:))]; %stim2 bounding box
box1 = [xyMatrix(1,2), xyMatrix(2,2), ...
    RectWidth(stimRect_shifted(1,:)), RectHeight(stimRect_shifted(1,:))]; %stim1 bounding box

overlaps = bboxOverlapRatio(box1, box2, 'Min');
for i = 2:length(weighted_x_trimmed2)
    stimRect_shifted2(i,:) = OffsetRect(stimRect_shifted2(i-1,:), weighted_x_trimmed2(i), weighted_y_trimmed2(i));
    [xyMatrix2(1,i+1),xyMatrix2(2,i+1)] = RectCenter(stimRect_shifted2(i,:)); % xyMatrix is 1 ahead of stimRect_shifted bc it includes starting pos
    
    box2 = [xyMatrix2(1,i+1), xyMatrix2(2,i+1), ...
        RectWidth(stimRect_shifted2(i,:)), RectHeight(stimRect_shifted2(i,:))]; %stim2 bounding box
    box1 = [xyMatrix(1,i+1), xyMatrix(2,i+1), ...
        RectWidth(stimRect_shifted(i,:)), RectHeight(stimRect_shifted(i,:))]; %stim1 bounding box

    if sum(stimRect_shifted(i,:) < 0)>=1 || sum(stimRect_shifted2(i,:) < 0)>=1
        overlaps(i) = 0;
    else
        overlaps(i) = bboxOverlapRatio(box1, box2, 'Min');
    end
end

%% keep updating if the stim1/stim2 overlap falls under .5
lessThanHalf = find(overlaps<.75);
startRegen = GetSecs();
while length(lessThanHalf)>=1 
    disp('Overlap fell below .75, regenerating...')
    newRands = 2*pi*rand([1,length(lessThanHalf)]); %one for stim2
    s2_newdel_x2 = speed*cos(2*pi*newRands); %only as many as we need
    s2_newdel_y2 = speed*sin(2*pi*newRands);
    
    newRands = 2*pi*rand([1,length(lessThanHalf)]); %one for stim1
    s1_newdel_x2 = speed*cos(2*pi*newRands); %only as many as we need
    s1_newdel_y2 = speed*sin(2*pi*newRands);
    
    for k = 1:length(lessThanHalf) %reinsert new random thetas into our original vector
        s2_delta_x2(lessThanHalf(k)) = s2_newdel_x2(k);
        s2_delta_y2(lessThanHalf(k)) = s2_newdel_y2(k);
    end
    
    for k = 1:length(lessThanHalf) %taken from gentraj1 for stim1
        delta_x(lessThanHalf(k)) = s1_newdel_x2(k);
        delta_y(lessThanHalf(k)) = s1_newdel_y2(k);
    end
    
    %stim2
    weighted_x2 = s2_delta_x2*weightMatrix2'; %transposed matrix to put oldest on top and help calculate in sequential order
    weighted_y2 = s2_delta_y2*weightMatrix2';

    weighted_x_trimmed2 = weighted_x2((length(weights)):end);
    weighted_y_trimmed2 = weighted_y2((length(weights)):end);
    
    %stim1
    weighted_x = delta_x*weightMatrix'; %transposed matrix to put oldest on top and help calculate in sequential order
    weighted_y = delta_y*weightMatrix';
    
    weighted_x_trimmed = weighted_x((length(weights)):end);
    weighted_y_trimmed = weighted_y((length(weights)):end);
    
    %%%
    stimRect_shifted2 = OffsetRect(stimRect, weighted_x2(1), weighted_y2(1));
    stimRect_shifted = OffsetRect(stimRect, weighted_x(1), weighted_y(1));
  
    [xyMatrix2(1,2),xyMatrix2(2,2)] = RectCenter(stimRect_shifted2(1,:)); % (1,:) = x, (2,:) = y
    [xyMatrix(1,2),xyMatrix(2,2)] = RectCenter(stimRect_shifted(1,:)); % (1,:) = x, (2,:) = y
    
    
    box2 = [xyMatrix2(1,2), xyMatrix2(2,2), ...
        RectWidth(stimRect_shifted2(1,:)), RectHeight(stimRect_shifted2(1,:))]; %stim2 bounding box
    box1 = [xyMatrix(1,2), xyMatrix(2,2), ...
        RectWidth(stimRect_shifted(1,:)), RectHeight(stimRect_shifted(1,:))]; %stim1 bounding box
    overlaps = bboxOverlapRatio(box1,box2, 'Min');
    
    for i = 2:length(weighted_x_trimmed2)
        %re-shift stim1
        stimRect_shifted(i,:) = OffsetRect(stimRect_shifted(i-1,:), weighted_x_trimmed(i), weighted_y_trimmed(i));
        [xyMatrix(1,i+1),xyMatrix(2,i+1)] = RectCenter(stimRect_shifted(i,:));
        %re-shift stim2
        stimRect_shifted2(i,:) = OffsetRect(stimRect_shifted2(i-1,:), weighted_x_trimmed2(i), weighted_y_trimmed2(i));
        [xyMatrix2(1,i+1),xyMatrix2(2,i+1)] = RectCenter(stimRect_shifted2(i,:));
        
        box2 = [xyMatrix2(1,i+1), xyMatrix2(2,i+1), ...
            RectWidth(stimRect_shifted2(i,:)), RectHeight(stimRect_shifted2(i,:))]; %stim2 bounding box
        box1 = [xyMatrix(1,i+1), xyMatrix(2,i+1), ...
            RectWidth(stimRect_shifted(i,:)), RectHeight(stimRect_shifted(i,:))]; %stim1 bounding box

        if sum(stimRect_shifted(i,:) < 0)>=1 || sum(stimRect_shifted2(i,:) < 0)>=1
            overlaps(i) = 0; % want to make sure that it doesn't go out of bounds
        else
            overlaps(i) = bboxOverlapRatio(box1, box2, 'Min');
        end
    end
    
    lessThanHalf = find(overlaps<.75);
end


%% convert back into weighted theta values

weighted_theta2 = atan(weighted_y_trimmed2./weighted_x_trimmed2);
weighted_theta = atan(weighted_y_trimmed./weighted_x_trimmed);