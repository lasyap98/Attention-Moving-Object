function validIdx = CollectValidTraj(possTrials, borderRect)

validIdx = [];

for j = 1:length(possTrials)
    start_xs = possTrials{j}(:,1);
    start_ys = possTrials{j}(:,2);
    end_xs = possTrials{j}(:,3);
    end_ys = possTrials{j}(:,4);
    
    if min(start_xs)>= borderRect(1) && min(start_ys)>= borderRect(2) ...
            && max(end_xs)<= borderRect(3) && max(end_ys)<= borderRect(4)
        validIdx(end+1) = j; %record the indices in possTrials where trajectory stays within border
    end
end

end