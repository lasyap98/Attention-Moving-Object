

% clear lag stim1_x stim1_y stim2_x stim2_y %currTrial_stim currTrial_stim2



%% organize all x and y values by trial
for b = 1:s %1 = face, 2 = house
for trial = 1:Lag.tperS
% eyetracking-x values across time
    lag.eye_x{b}{trial} = Lag.EyeData{b}.mx{trial}; %get values after fixating at center for 500ms%
end
for trial = 1:Lag.tperS
% eyetracking-y values across time
    lag.eye_y{b}{trial} = Lag.EyeData{b}.my{trial}; %get values after fixating at center for 500ms%
end
for trial = 1:Lag.tperS
% eyetracking-a values across time
    lag.eye_a{b}{trial} = Lag.EyeData{b}.ma{trial}; %get values after fixating at center for 500ms%
end
end

% stimulus1-x values across time
for trial = 1:length(Lag.currTrial_stim)
    lag.face_x{trial} = [];
    for i = 1:length(Lag.currTrial_stim{trial})
        lag.face_x{trial}(end+1)= mean([Lag.currTrial_stim{trial}(i,1),Lag.currTrial_stim{trial}(i,3)]);
    end
end

% stimulus1-y values across time
for trial = 1:length(Lag.currTrial_stim)
    lag.face_y{trial} = [];
    for i = 1:length(Lag.currTrial_stim{trial})
        lag.face_y{trial}(end+1)= mean([Lag.currTrial_stim{trial}(i,2),Lag.currTrial_stim{trial}(i,4)]);
    end
end

% stimulus2-x values across time
for trial = 1:length(Lag.currTrial_stim2)
    lag.house_x{trial} = [];
    for i = 1:length(Lag.currTrial_stim2{trial})
        lag.house_x{trial}(end+1)= mean([Lag.currTrial_stim2{trial}(i,1),Lag.currTrial_stim2{trial}(i,3)]);
    end
end

% stimulus2-y values across time
for trial = 1:length(Lag.currTrial_stim2)
    lag.house_y{trial} = [];
    for i = 1:length(Lag.currTrial_stim2{trial})
        lag.house_y{trial}(end+1)= mean([Lag.currTrial_stim2{trial}(i,2),Lag.currTrial_stim2{trial}(i,4)]);
    end
end

%% remove frames that had no eye data
for b = 1:s
for trial = 1:Lag.tperS
    zero_pupil1= find(lag.eye_a{b}{trial} ==0);
    
    beforeFix = find(Lag.EyeData{b}.FixDoneT{trial}<Lag.flipTimes{b}{trial}(2)); 

    if lag.eye_a{b}{trial}(beforeFix(end-1))==0
        %for x-y we don't want to remove the same ones as eye-a in case there was 
        zero_pupilb4fix = intersect(zero_pupil1, beforeFix); %find when zero pupil before fixation 
        nozero_b4fix = beforeFix(1:(end-length(zero_pupilb4fix)+1)); %remove before-fixation points to account for no pupil frames, use for eye values
    else
        zero_pupilb4fix = intersect(zero_pupil1, beforeFix); %find when zero pupil before fixation 
        nozero_b4fix = beforeFix(1:(end-length(zero_pupilb4fix))); %remove before-fixation points to account for no pupil frames, use for eye values
    end
    
    lag.eye_x{b}{trial} = lag.eye_x{b}{trial}(nozero_b4fix(end-1):end);
    lag.eye_y{b}{trial} = lag.eye_y{b}{trial}(nozero_b4fix(end-1):end);
    lag.eye_a{b}{trial} = lag.eye_a{b}{trial}(beforeFix(end-1):end);
        
    %where do we not see the pupil? (after fixation now)
    zero_pupil2 = find(lag.eye_a{b}{trial}==0);
        
    del = 0;
    
    if b==1 %face
        for frame = zero_pupil2
            lag.face_x{trial}(frame-del) = [];
            lag.face_y{trial}(frame-del) = [];
            del = del+1;
        end
        assert(length(lag.face_x{trial})==length(lag.eye_x{b}{trial}), 'Array sizes don''t match! (face)')
    else
        for frame = zero_pupil2
            lag.house_x{trial}(frame-del) = [];
            lag.house_y{trial}(frame-del) = [];
            del = del+1;
        end
        
        assert(length(lag.house_x{trial})==length(lag.eye_x{b}{trial}), 'Array sizes don''t match! (house)')
    end
   
end
end

%% find max correlation for optimal lag

for b = 1:s
    
for t = 1:Lag.tperS
    
    for n=1:30
        
        s1x_temp = lag.face_x{t}(1:end+1-n);
        s1y_temp = lag.face_y{t}(1:end+1-n);
        
        s2x_temp = lag.house_x{t}(1:end+1-n);
        s2y_temp = lag.house_y{t}(1:end+1-n);

        ex_temp = lag.eye_x{b}{t}(n:end);
        ey_temp = lag.eye_y{b}{t}(n:end);
        
        if b==1
            r_score.face{t}(n,1) = corr2(s1x_temp, ex_temp);
            r_score.face{t}(n,2) = corr2(s1y_temp, ey_temp);
        else
            r_score.house{t}(n,1) = corr2(s2x_temp, ex_temp);
            r_score.house{t}(n,2) = corr2(s2y_temp, ey_temp);
        end
    end
    
    if b==1
        max_corr.face(t,1) = max(r_score.face{t}(:,1)); % max correlation for x values
        max_corr.face(t,2) = max(r_score.face{t}(:,2)); % max correlation for y values
        
        max_lag.face(t,1) = find(r_score.face{t}(:,1)==max(r_score.face{t}(:,1))); % optimal lag for x
        max_lag.face(t,2) = find(r_score.face{t}(:,2)==max(r_score.face{t}(:,2))); % optimal lag for y
    else
        max_corr.house(t,1) = max(r_score.house{t}(:,1)); % max correlation for x values
        max_corr.house(t,2) = max(r_score.house{t}(:,2)); % max correlation for y values
        
        max_lag.house(t,1) = find(r_score.house{t}(:,1)==max(r_score.house{t}(:,1))); % optimal lag for x
        max_lag.house(t,2) = find(r_score.house{t}(:,2)==max(r_score.house{t}(:,2))); % optimal lag for y
    end
end


end

%% plot correlation scores by lag
% 
% iTrial = 1;
% 
% subplot(1,2,1)
% plot(r_score.face{iTrial}(:,1))
% hold on
% plot(r_score.face{iTrial}(:,2))
% legend('Stim 1 (x)', 'Stim 1 (y)')
%     
% subplot(1,2,2)
% plot(r_score.house{iTrial}(:,1))
% hold on
% plot(r_score.house{iTrial}(:,2))
% legend('Stim 2 (x)', 'Stim 2 (y)')

%% collect average lag across the lag trials
% to use in regular trials

lag.avgmaxlag.face_x = round(mean(max_lag.face(:,1)));
lag.avgmaxlag.face_y = round(mean(max_lag.face(:,2)));
lag.avgmaxlag.house_x = round(mean(max_lag.house(:,1)));
lag.avgmaxlag.house_y = round(mean(max_lag.house(:,2)));

lag.avgmaxlag.face = round(mean(mean(max_lag.face())));
lag.avgmaxlag.house = round(mean(mean(max_lag.house())));
lag.avgmaxlag.facehouse = round(mean([lag.avgmaxlag.face, lag.avgmaxlag.house]));

% save(sprintf('analyses_%s.mat', Params.sub.info), 'lag')
%% sliding window of correlation (within lag trials)

% win_size = 1; %sec
% win_frames = round(win_size/refresh);
% 
% % EACH trial and EACH stim has an optimal lag
% for t = 1:Lag.tperS
%     
% %     %1-1 lineup of eye and stim-1 and stim-2
% %     x1_trial = lag.eye_x{1}{t}(max_lag.face(t,1):end);
% %     s1x_trial = lag.face_x{t}(1:end-max_lag.face(t,1)+1);
% %     y1_trial = lag.eye_y{1}{t}(max_lag.face(t,2):end);
% %     s1y_trial = lag.face_y{t}(1:end-max_lag.face(t,2)+1);
%     
%     x2_trial = lag.eye_x{2}{t}(max_lag.house(t,1):end);
%     s2x_trial = lag.house_x{t}(1:end-max_lag.house(t,1)+1);
%     y2_trial = lag.eye_y{2}{t}(max_lag.house(t,2):end);
%     s2y_trial = lag.house_y{t}(1:end-max_lag.house(t,2)+1);
%     
% %     win_shifts_s1x = length(s1x_trial)-win_frames+1;
%     win_shifts_s2x = length(s2x_trial)-win_frames+1;
% %     win_shifts_s1y = length(s1y_trial)-win_frames+1;
%     win_shifts_s2y = length(s2y_trial)-win_frames+1;
%   
% %     lag.corr_win.facex{t} = zeros([1, win_shifts_s1x]);
% %     lag.corr_win.facey{t} = zeros([1, win_shifts_s1y]);
%     lag.corr_win.housex{t} = zeros([1, win_shifts_s2x]);
%     lag.corr_win.housey{t} = zeros([1, win_shifts_s2y]);
%     
% %     for w = 1:win_shifts_s1x        
% %         lag.corr_win.facex{t}(w) = corr2(x1_trial(w:w+win_frames-1), s1x_trial(w:w+win_frames-1));
% %     end
%     for w = 1:win_shifts_s2x        
%         lag.corr_win.housex{t}(w) = corr2(x2_trial(w:w+win_frames-1), s2x_trial(w:w+win_frames-1));
%     end
% %     for w = 1:win_shifts_s1y        
% %         lag.corr_win.facey{t}(w) = corr2(y1_trial(w:w+win_frames-1), s1y_trial(w:w+win_frames-1));
% %     end
%     for w = 1:win_shifts_s2y        
%         lag.corr_win.housey{t}(w) = corr2(y2_trial(w:w+win_frames-1), s2y_trial(w:w+win_frames-1));
%     end
%     
%     
% end


%%
% % NOW GO BACK TO MAIN TRIAL DATA (one block only)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if subject==3
%     block = [2, 0, 1, 2, 0, 1];
% elseif subject==1
%     block = [1, 2, 0, 1, 2, 0]; %face, house
% elseif subject==2
%     block = [0, 1, 2, 0, 1, 2]; %both, face, house
% end
% 
% face.eyedata = Block_EyeData(block==1);
% house.eyedata = Block_EyeData(1);
% both.eyedata = Block_EyeData(block==0);
% 
% face.t2Run = Block_t2Run(block==1);
% house.t2Run = Block_t2Run(1);
% both.t2Run = Block_t2Run(block==0);
% 
% for b = 1:2
%     face.face_stim{b} = pT(face.t2Run{b});
%     face.house_stim{b} = pT2(face.t2Run{b});
%     
%     house.face_stim{b} = pT(house.t2Run{b});
%     house.house_stim{b} = pT2(house.t2Run{b});
%     
%     both.face_stim{b} = pT(both.t2Run{b});
%     both.house_stim{b} = pT2(both.t2Run{b});
% end
% 
% %% organize all x and y values by trial - face
% 
% for b = 1:2
%     
% for trial = 1:trialNum
% % eyetracking-x, y, a values across time
% %     face.eye_x{b}{trial} = face.eyedata{b}.mx{trial}; %30+ get values after fixating at center for 500ms%
% %     face.eye_y{b}{trial} = face.eyedata{b}.my{trial}; %30+ get values after fixating at center for 500ms%
% %     face.eye_a{b}{trial} = face.eyedata{b}.ma{trial}; %30+ get values after fixating at center for 500ms%
% %     
%     house.eye_x{b}{trial} = house.eyedata{b}.mx{trial}; %30+ get values after fixating at center for 500ms%
%     house.eye_y{b}{trial} = house.eyedata{b}.my{trial}; %30+ get values after fixating at center for 500ms%
%     house.eye_a{b}{trial} = house.eyedata{b}.ma{trial}; %30+ get values after fixating at center for 500ms%
%     
% %     both.eye_x{b}{trial} = both.eyedata{b}.mx{trial}; %30+ get values after fixating at center for 500ms%
% %     both.eye_y{b}{trial} = both.eyedata{b}.my{trial}; %30+ get values after fixating at center for 500ms%
% %     both.eye_a{b}{trial} = both.eyedata{b}.ma{trial}; %30+ get values after fixating at center for 500ms%
% end
% 
% % face x,y values across time
% for trial = 1:length(face.face_stim{b})
%     face.face_x{b}{trial} = [];
%     face.face_y{b}{trial} = [];
%     face.house_x{b}{trial} = [];
%     face.house_y{b}{trial} = [];
%     for i = 1:length(face.face_stim{b}{trial})
%         face.face_x{b}{trial}(end+1)= mean([face.face_stim{b}{trial}(i,1),face.face_stim{b}{trial}(i,3)]);
%         face.face_y{b}{trial}(end+1)= mean([face.face_stim{b}{trial}(i,2),face.face_stim{b}{trial}(i,4)]);
%         face.house_x{b}{trial}(end+1)= mean([face.house_stim{b}{trial}(i,1),face.house_stim{b}{trial}(i,3)]);
%         face.house_y{b}{trial}(end+1)= mean([face.house_stim{b}{trial}(i,2),face.house_stim{b}{trial}(i,4)]);
%     end
% end
% 
% for trial = 1:length(house.face_stim{b})
%     house.face_x{b}{trial} = [];
%     house.face_y{b}{trial} = [];
%     house.house_x{b}{trial} = [];
%     house.house_y{b}{trial} = [];
%     for i = 1:length(house.face_stim{b}{trial})
%         house.face_x{b}{trial}(end+1)= mean([house.face_stim{b}{trial}(i,1),house.face_stim{b}{trial}(i,3)]);
%         house.face_y{b}{trial}(end+1)= mean([house.face_stim{b}{trial}(i,2),house.face_stim{b}{trial}(i,4)]);
%         house.house_x{b}{trial}(end+1)= mean([house.house_stim{b}{trial}(i,1),house.house_stim{b}{trial}(i,3)]);
%         house.house_y{b}{trial}(end+1)= mean([house.house_stim{b}{trial}(i,2),house.house_stim{b}{trial}(i,4)]);
%     end
% end
% 
% for trial = 1:length(both.face_stim{b})
%     both.face_x{b}{trial} = [];
%     both.face_y{b}{trial} = [];
%     both.house_x{b}{trial} = [];
%     both.house_y{b}{trial} = [];
%     for i = 1:length(both.face_stim{b}{trial})
%         both.face_x{b}{trial}(end+1)= mean([both.face_stim{b}{trial}(i,1),both.face_stim{b}{trial}(i,3)]);
%         both.face_y{b}{trial}(end+1)= mean([both.face_stim{b}{trial}(i,2),both.face_stim{b}{trial}(i,4)]);
%         both.house_x{b}{trial}(end+1)= mean([both.house_stim{b}{trial}(i,1),both.house_stim{b}{trial}(i,3)]);
%         both.house_y{b}{trial}(end+1)= mean([both.house_stim{b}{trial}(i,2),both.house_stim{b}{trial}(i,4)]);
%     end
% end
% end
% %% remove frames that had no eye data
% 
% for b = 1:2
%     
% for trial = 1:trialNum
%     zero_pupil1= find(house.eye_a{b}{trial} ==0);
%     
%     beforeFix = find(EyeData.FixDoneT{trial}<flipTimes{trial}(2)); 
% 
%     if house.eye_a{b}{trial}(beforeFix(end-1))==0
%         %for x-y we don't want to remove the same ones as eye-a in case there was 
%         zero_pupilb4fix = intersect(zero_pupil1, beforeFix); %find when zero pupil before fixation 
%         nozero_b4fix = beforeFix(1:(end-length(zero_pupilb4fix)+1)); %remove before-fixation points to account for no pupil frames, use for eye values
%     else
%         zero_pupilb4fix = intersect(zero_pupil1, beforeFix); %find when zero pupil before fixation 
%         nozero_b4fix = beforeFix(1:(end-length(zero_pupilb4fix))); %remove before-fixation points to account for no pupil frames, use for eye values
%     end
%     
%     house.eye_x{b}{trial} = house.eye_x{b}{trial}(nozero_b4fix(end-1):end);
%     house.eye_y{b}{trial} = house.eye_y{b}{trial}(nozero_b4fix(end-1):end);
%     house.eye_a{b}{trial} = house.eye_a{b}{trial}(beforeFix(end-1):end);
%         
%     %where do we not see the pupil? (after fixation now)
%     zero_pupil2 = find(house.eye_a{b}{trial}==0);
%         
%     del = 0;
%     
%     for frame = zero_pupil2
%         house.face_x{b}{trial}(frame-del) = [];
%         house.face_y{b}{trial}(frame-del) = [];
%         house.house_x{b}{trial}(frame-del) = [];
%         house.house_y{b}{trial}(frame-del) = [];
%         
%         del = del+1;
%     end
% 
%     assert(length(house.house_x{b}{trial})==length(house.eye_x{b}{trial}), 'Array sizes don''t match!')
%     
% end
% end
% %% sliding window of correlation (previously known as version #2)
% 
% win_size = 1; %sec
% win_frames = round(win_size/refresh);
% 
% % EACH trial and EACH stim has an optimal lag
% 
% for t = 1:trialNum
%     
%     %1-1 lineup of eye-x and stim-1 and stim-2
%     x1_trial = house.eye_x{b}{t}(lag.avgmaxlag.house:end);
%     s1x_trial = house.face_x{b}{t}(1:end-lag.avgmaxlag.house+1);
%     x2_trial = house.eye_x{b}{t}(lag.avgmaxlag.house:end);
%     s2x_trial = house.house_x{b}{t}(1:end-lag.avgmaxlag.house+1);
% 
%     y1_trial = house.eye_y{b}{t}(lag.avgmaxlag.house:end);
%     s1y_trial = house.face_y{b}{t}(1:end-lag.avgmaxlag.house+1);
%     y2_trial = house.eye_y{b}{t}(lag.avgmaxlag.house:end);
%     s2y_trial = house.house_y{b}{t}(1:end-lag.avgmaxlag.house+1);
% 
%     win_shifts_s1x = length(s1x_trial)-win_frames+1;
%     win_shifts_s1y = length(s1y_trial)-win_frames+1;
%     win_shifts_s2x = length(s2x_trial)-win_frames+1;
%     win_shifts_s2y = length(s2y_trial)-win_frames+1;
%   
%     corr_win.stim1x{t} = zeros([1, win_shifts_s1x]);
%     corr_win.stim1y{t} = zeros([1, win_shifts_s1y]);
%     corr_win.stim2x{t} = zeros([1, win_shifts_s2x]);
%     corr_win.stim2y{t} = zeros([1, win_shifts_s2y]);
%     
%     for w = 1:win_shifts_s1x        
%         corr_win.stim1x{t}(w) = corr2(x1_trial(w:w+win_frames-1), s1x_trial(w:w+win_frames-1));
%     end
%     for w = 1:win_shifts_s1y        
%         corr_win.stim1y{t}(w) = corr2(y1_trial(w:w+win_frames-1), s1y_trial(w:w+win_frames-1));
%     end
%     for w = 1:win_shifts_s2x        
%         corr_win.stim2x{t}(w) = corr2(x2_trial(w:w+win_frames-1), s2x_trial(w:w+win_frames-1));
%     end
%     for w = 1:win_shifts_s2y        
%         corr_win.stim2y{t}(w) = corr2(y2_trial(w:w+win_frames-1), s2y_trial(w:w+win_frames-1));
%     end 
%     
%     corrAnalyses.x1{t} = x1_trial; %eye-x values shifted by optimal lag for stim1
%     corrAnalyses.y1{t} = y1_trial; %eye-y values shifted by optimal lag for stim1
% 
%     corrAnalyses.x2{t} = x2_trial; %eye-x values shifted by optimal lag for stim2
%     corrAnalyses.y2{t} = y2_trial; %eye-y values shifted by optimal lag for stim2
%     
%     corrAnalyses.s1x{t} = s1x_trial;
%     corrAnalyses.s1y{t} = s1y_trial;
% 
%     corrAnalyses.s2x{t} = s2x_trial;
%     corrAnalyses.s2y{t} = s2y_trial;
% 
%     mean_face(1,t) = mean(corr_win.stim1x{t});
%     mean_face(2,t) = mean(corr_win.stim1y{t});
% 
%     mean_house(1,t) = mean(corr_win.stim2x{t});
%     mean_house(2,t) = mean(corr_win.stim2y{t});
% end
%     
% mean(mean_face(1,:))
% mean(mean_face(2,:))
% mean(mean_house(1,:))
% mean(mean_house(2,:))
% 
% %% temporal correlation: stim1 v stim 2
% 
% win_size = 1; %sec
% win_frames = round(win_size/refresh);
% 
% for t = 1:length(stim1_x)
% 
%     win_shifts_s1s2 = length(stim2_x{t})-win_frames+1; %stim1 and stim2 will have same length, x or y
%     
%     corr_win.s1s2x{t} = zeros([1, win_shifts_s1s2]);
%     corr_win.s1s2y{t} = zeros([1, win_shifts_s1s2]);
%     
%     for w = 1:win_shifts_s1s2        
%         corr_win.s1s2x{t}(w) = corr2(stim1_x{t}(w:w+win_frames-1), stim2_x{t}(w:w+win_frames-1));
%     end 
%     for w = 1:win_shifts_s1s2        
%         corr_win.s1s2y{t}(w) = corr2(stim1_y{t}(w:w+win_frames-1), stim2_y{t}(w:w+win_frames-1));
%     end  
% 
% end
% 
% mean_s1s2corr = zeros([2, trialNum]);
% 
% for t = 1:length(corr_win.s1s2x)
%     mean_s1s2corr(1,t) = mean(corr_win.s1s2x{t});
%     mean_s1s2corr(2,t) = mean(corr_win.s1s2y{t});
% end

%%
% 
% 
% mean(mean_s1s2corr(1,:))
% mean(mean_s1s2corr(2,:))
% 
% %% plot correlation window analysis #2
% 
% % subplot(2,2,1)
% % plot(corrAnalyses.x1{iTrial});
% % hold on
% % plot(corrAnalyses.s1x{iTrial});
% % hold on
% % plot(corrAnalyses.s2x{iTrial});
% % legend('Gaze(x)','Attended Stimulus(x)','Unattended Stimulus(x)');
% % title(sprintf("Trial: %d", iTrial));
% % % %
% % subplot(2,2,2)
% iTrial = 1;
% plot(corr_win.s1s2x{iTrial});
% hold on
% plot(corr_win.stim1x{iTrial}, 'Color',[0.8500, 0.3250, 0.0980]);
% hold on
% plot(corr_win.stim2x{iTrial},'Color',[0.9290, 0.6940, 0.1250]);
% legend('S1 v S2 Corr(x)', 'Stim1 Corr(x)','Stim2 Corr(x)');
% title(sprintf("Trial: %d, Window size: %d second(s)", iTrial, win_size));
% 
% % % % % % % % %
% iTrial = idxR(2);
% subplot(2,2,3)
% plot(corrAnalyses2.x1{iTrial});
% hold on
% plot(corrAnalyses2.s1x{iTrial});
% hold on
% plot(corrAnalyses2.s2x{iTrial});
% legend('Gaze(x)','Attended Stimulus(x)','Unattended Stimulus(x)');
% title(sprintf("Trial: %d", iTrial));
% %
% subplot(2,2,4)
% plot(corr_win2.stim1.x{1, iTrial}, 'Color',[0.8500, 0.3250, 0.0980]);
% hold on
% plot(corr_win2.stim2.x{1, iTrial},'Color',[0.9290, 0.6940, 0.1250]);
% legend('Stim1 Corr(x)','Stim2 Corr(x)');
% title(sprintf("Trial: %d, Window size: %d second(s)", iTrial, win_size));
% 
% subplot(2,2,3)
% plot(corrAnalyses2.y1{iTrial});
% hold on
% plot(corrAnalyses2.s1y{iTrial});
% hold on
% plot(corrAnalyses2.s2y{iTrial});
% legend('Gaze(y)','Attended Stimulus(y)','Unattended Stimulus(y)');
% title(sprintf("Trial: %d", iTrial));
% %
% subplot(2,2,4)
% plot(corr_win.s1s2.y{iTrial});
% hold on
% plot(corr_win.stim1.y{1, iTrial}, 'Color',[0.8500, 0.3250, 0.0980]);
% hold on
% plot(corr_win.stim2.y{1, iTrial},'Color',[0.9290, 0.6940, 0.1250]);
% legend('S1 v S2 Corr(y)', 'Stim1 Corr(y)','Stim2 Corr(y)');
% title(sprintf("Trial: %d, Window size: %d second(s)", iTrial, win_size));
% 
