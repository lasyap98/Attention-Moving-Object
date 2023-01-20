global Params Display
%% Generates random trajectory for ONE trial

x = linspace(-3,3,64);
sd = 1;
y = exp(-0.5*(x/sd).^2);
weights = (exp(-0.5*(x/sd).^2))/4; % to be assigned to 'n=length' previous directions
%% mean theta will need to be weighted average of n previous directions
% collect past n directions to include in weighted average
% 'newest' direction on top (1)

speed = 1; % in pixels: the bigger the stimulus, the faster it will need to go   
total_frames = round(Params.trialLength/Display.refresh); % will be based on refresh rate
randn = rand([1,total_frames+length(weights)-1]); %how many random thetas do we want
randTheta = 2*pi*randn; %create the random thetas

%% create matrix of weights, should have full table adjusted for each frame
% can use for loop
% the greatest weight should be assigned to the most recent frame and shift accordingly

weightMatrix = zeros(length(randTheta));

for i = 1:length(randTheta)
    for k = 0:length(weights)-1
        if (i+k) > length(randTheta)
            break
        end
        weightMatrix(i+k,i) = weights(k+1); 
    end
end        
    
%% convert thetas to delta-x and delta-y values

delta_x = speed*cos(randTheta);
delta_y = speed*sin(randTheta);
    
%% multiply matrices to get weighted delta-x and delta-y values
% 
weighted_x = delta_x*weightMatrix'; %transposed matrix to put oldest on top and help calculate in sequential order
weighted_y = delta_y*weightMatrix';

% max_trans = 2; %restrict translation to a max of n pixels
% translation = sqrt(weighted_x.^2+weighted_y.^2);
%translation(translation>max_trans) = max_trans; %where does translation exceed

% while sum(translation>max_trans)
%     randn = rand([1,sum(translation>max_trans)]); %how many random thetas do we want
%     randTheta = 2*pi*randn;
%     
%     for k = 1:length(find(translation>max_trans))
%         delta_x(k) = speed*cos(randTheta(k));
%         delta_y(k) = speed*sin(randTheta(k));
%     end
%     
%     weighted_x = delta_x*weightMatrix'; %transposed matrix to put oldest on top and help calculate in sequential order
%     weighted_y = delta_y*weightMatrix';
%     
%     translation = sqrt(weighted_x.^2+weighted_y.^2);
% 
% end


% only keep the ones with full weights calculated
weighted_x_trimmed = weighted_x((length(weights)):end);
weighted_y_trimmed = weighted_y((length(weights)):end);

%% update position matrix

xyMatrix = zeros([2,length(weighted_x_trimmed)]); % (1,:) = x, (2,:) = y, 
% xyMatrix will have one more than stimRect_shifted bc we record starting position (xy center)
xyMatrix(:,1) = [stim_xCtr;stim_yCtr]; % starting position
stimRect_shifted = OffsetRect(stimRect, weighted_x(1), weighted_y(1));
[xyMatrix(1,2),xyMatrix(2,2)] = RectCenter(stimRect_shifted(1,:));
for i = 2:length(weighted_x_trimmed)
    stimRect_shifted(i,:) = OffsetRect(stimRect_shifted(i-1,:), weighted_x_trimmed(i), weighted_y_trimmed(i));
    [xyMatrix(1,i+1),xyMatrix(2,i+1)] = RectCenter(stimRect_shifted(i,:));
end


%% convert back into weighted theta values

weighted_theta = atan(weighted_y./weighted_x);