global el EyeData currPosID Display lagOn

if ET
    errormes=Eyelink('CheckRecording');
%         if(errormes~=0)
%         end
    if Eyelink('NewFloatSampleAvailable') > 0
        % get the sample in the form of an event structure
        evt = Eyelink('NewestFloatSample');
        x = evt.gx(evt.gx~=el.MISSING_DATA);
        y = evt.gy(evt.gy~=el.MISSING_DATA);
        a = evt.pa(evt.pa~=-el.MISSING_DATA);
        tm = evt.time; %sync with eyelink time events
        fixT = GetSecs(); %doesn't account for blinks, so fix this to get accurate timestamps
        
        if lagOn
            Lag.EyeData{s}.mx{iTrial}=[Lag.EyeData{s}.mx{iTrial} x];
            Lag.EyeData{s}.my{iTrial}=[Lag.EyeData{s}.my{iTrial} y];
            Lag.EyeData{s}.ma{iTrial}=[Lag.EyeData{s}.ma{iTrial} a];
            Lag.EyeData{s}.time{iTrial}=[Lag.EyeData{s}.time{iTrial} tm];
            Lag.EyeData{s}.FixDoneT{iTrial} = [Lag.EyeData{s}.FixDoneT{iTrial} fixT]; 
        else
            EyeData.mx{iTrial}=[EyeData.mx{iTrial} x];
            EyeData.my{iTrial}=[EyeData.my{iTrial} y];
            EyeData.ma{iTrial}=[EyeData.ma{iTrial} a];
            EyeData.time{iTrial}=[EyeData.time{iTrial} tm];
            EyeData.FixDoneT{iTrial} = [EyeData.FixDoneT{iTrial} fixT]; % record each time we get valid eyetracking data point
        end
        % check whether fix the center
        if ~ (isempty(x) || isempty(y) || isempty(a))
            if x~=el.MISSING_DATA && y~=el.MISSING_DATA && a>0 % do we have valid data and is the pupil visible?
                dist_from_fix = round(( (x-(Display.rect(3)/2))^2 + (y-(Display.rect(4)/2))^2 ) ^ (1/2));

                if dist_from_fix <= 50
                    Fixated = 1;
                    
                    if currPosID==1 % until the motion turns on, keep noting the difference in time from inital start to 'now'
                        if lagOn
                            Lag_InitialFixDur(iTrial) = fixT-Lag_EyeStart(iTrial);
                        else
                            InitialFixDur(iTrial) = fixT-EyeStart(iTrial);
                        end
                    end
                else
                    Fixated = 0;
                end
            end
        else
            Fixated = 0;
        end
    end
else
    fixBreak=0;
end
