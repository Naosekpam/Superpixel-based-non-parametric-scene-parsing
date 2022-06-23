
RunSystemMake3D;
stage = 1;
GT = cell(length(testList),1);
Y = cell(length(testList),1);
for i = 1:length(testList)
    
    i
    query = testList{i};
    load([outputDir query  '.mat']);%num2str(stage)
    load([labelDir query '.mat']); %Assume that the GT of the test file is in the TrainingLabel di
    valid = S~=0; 
    GT{i} = S(valid);
    Y{i} = predictLabel(valid);
    
end

GT = cell2mat(GT);
Y = cell2mat(Y);

pixelAcc = sum(GT==Y)/length(GT)

accuracies = zeros(NC,1);
for cc = 1 : NC
    if sum(GT==cc)~=0
        accuracies(cc) = sum((GT==Y)&(GT==cc))/sum(GT==cc);
    end
end
avAcc = mean(accuracies)
