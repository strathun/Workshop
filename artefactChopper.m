function [sortedData,breakStack,breakStackRev] = artefactChopper(data,Fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%LOG: have event detection working. Check stacker and then figure out how
%to partition into different channels!!

gap = ceil(Fs/68.2);
searchDistance = ceil(Fs/50);
meanData = mean(data);

TF1 = islocalmin(meanData, 'MinSeparation', searchDistance);
TF2 = find(TF1>0);
meanData2 = meanData(TF2);
minuse = min(meanData2);
meanData3 = find(meanData2<(minuse*6.5));
% eventVals = meanData2(meanData3);
eventsNew = TF2(meanData3);

chop = 1;
[channels, ~] = size(data);
[divs] = length(eventsNew);
sortedCell = cell(1,channels);
prevBreak = cell(channels,1);

for i = 1:channels
    start = 1;
    for ii = 1:divs
        stop = eventsNew(ii);
        breakStackTemp = data(i,start:(stop-gap));
%         figure(chop)
%         plot(breakStackTemp)
        breakStack{chop} = {breakStackTemp};
        chop = chop + 1;
        start = stop + gap; 
    end
end


ii = 1;
for i = 1:length(breakStack)
    meanArray(ii,i-((ii-1)*divs)) = mean(cell2mat(breakStack{i}));
    stdArray(ii,i-((ii-1)*divs)) = std(cell2mat(breakStack{i}));
    sumArray(ii,i-((ii-1)*divs)) = sum(cell2mat(breakStack{i}));
    breakStackRev{ii,i-((ii-1)*divs)} = breakStack{i};
    if mod(i , divs) == 0   %breaks up Arrays into channels to simply sorting below
        ii = ii + 1;
    end
end


for i = 1:channels
    sort1{1,i} = meanArray(i,1);
    sort1{2,i} = [];
end

totaldivs = length(breakStack);

for i = 1:totaldivs
    
    %comparing current value to each channel
    for ii = 1:channels
        tempDiff(ii) = abs(meanArray(i) - sort1{1,ii});
    end
%     %comparing current value to each channel
%     for ii = 1:channels
%         if i > channels
%             tempDiff(ii) = abs(meanArray(i) - mean(prevBreak{ii}));
% %         tempDiff(ii) = abs(meanArray(i) - sort1{1,ii});
%         else
%             tempDiff(ii) = abs(meanArray(i) - sort1{1,ii});
%         end
%     end
    
    channelAssignment = min(tempDiff);
    channelAssignment = find(tempDiff == channelAssignment);
    sort1{2,channelAssignment} = [sort1{2,channelAssignment},meanArray(i)];
    sort1{1,channelAssignment} = mean(sort1{2,channelAssignment});   %updates channel rep to improve sorting accuracy
    sortedCell{1,channelAssignment} = [sortedCell{1,channelAssignment},breakStackRev{i}];
    prevBreak{channelAssignment} = cell2mat(breakStackRev{i});
end

for i = 1:channels
    [~, sortSizes(i)] = size(sortedCell{1,i});
end

sortedTemp = [];
for ii = 1:channels
    cellTemp = sortedCell{1,ii};
    for i = 1:sortSizes(ii)
        sortedTemp = [sortedTemp,cellTemp{i}];
    end
    sortTempHolder{ii} = sortedTemp;
    sortedTemp = [];
end

for i = 1:channels
    [~, sizerArray(i)] = size(sortTempHolder{i});
end

temper1 = min(sizerArray);
temper2 = find(sizerArray == temper1);

for i = 1:channels
    sortedDataTemp = sortTempHolder{i};
    sortedData(i,:) = sortedDataTemp(1:temper1);
end

end

