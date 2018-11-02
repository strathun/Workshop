function [sortedData,breakStack,breakStackRev,eventsCorrected] = artefactChopper_manual(data,Fs,ManAut)
%UNTITLED Summary of this function goes here
%   Inputs :
%       data : "raw" data. Should already be divided into number of
%               channels. These will likely not be correctly sorted.
%       Fs : Sampling rate for each individual channel.
%       ManAut : Type 'Manual' for manual sorting (recommended).
%   Outputs :
%       sortedData : raw data correctly sorted with reset artefacts
%                   removed.
%       breakStack :
%       breakStackRev : 
%       eventsCorrected : reset event times for reference outside of
%                         function. Indices are corrected for shifting of data after artefact
%                         removal. 
%%%%%%%%

%Starting definitions
gap = ceil(Fs/68.2);
searchDistance = ceil(Fs/50);
scrollGap = ceil(Fs/2);
meanData = mean(data);
[~, dataLength] = size(meanData) ; 
scrollTotal = floor(dataLength/scrollGap);
xdata = (1:1:dataLength) ;

    % Scrolls through data and lets user select resets on graph. Stores x -
    % coordinates as events.
    if ManAut == 'Manual'
        disp('Click all visible reset artefacts one time. Then press enter!')
        i = 1;
        eventsNow{1} = {} ; 
        while (i*scrollGap) <= dataLength
            xdataCurrent = xdata((( i - 1 ) * scrollGap )+1:( i * scrollGap ) ) ;
            ydataCurrent = meanData((( i - 1 ) * scrollGap )+1:( i * scrollGap ) ) ;
            plot(xdataCurrent, ydataCurrent)
            title([ num2str(i) ' of ' num2str(scrollTotal) ])
            [x, ~] = ginput;
            close
            eventsNow{1} = [eventsNow{1},x.'] ; 
            i = i+1 ; 
        end
    %     xdataCurrent = xdata((( i  ) * scrollGap )+1: end ) ;
    %     ydataCurrent = meanData((( i ) * scrollGap )+1: end ) ;
    %     plot(xdataCurrent, ydataCurrent)
    %     [x, ~] = ginput;
    %     eventsNow{1} = [eventsNow{1},x.'] ; 
        eventsNew = cell2mat(eventsNow{1});
        eventsNew = ceil(eventsNew);
            
    else
        %   Alternative automatic reset detection. 
        TF1 = islocalmin(meanData, 'MinSeparation', searchDistance);
        TF2 = find(TF1>0);
        meanData2 = meanData(TF2);
        minuse = min(meanData2);
        meanData3 = find(meanData2<(minuse*6.5));
        % eventVals = meanData2(meanData3);
        eventsNew = TF2(meanData3);
    end


%Pulls out sections between resets and stores them for sorting
    chop = 1;
    [channels, ~] = size(data);
    [divs] = length(eventsNew);
    sortedCell = cell(1,channels);
    prevBreak = cell(channels,1);  
    if channels > 4
    subTrace = 4;
    else
        subTrace = channels;
    end
    
    for i = 1:channels
        start = 1;
        for ii = 1:divs
            stop = eventsNew(ii);
            breakStackTemp = data(i,start:(stop-gap));
            breakStack{chop} = {breakStackTemp};
            chop = chop + 1;
            start = stop + gap; 
        end
    end


%corrects event index for gap deletions. Only relevant to exported data.
    for i = 1:length(eventsNew)
        eventsCorrected(i) = eventsNew(i) - ( ( ( i-1) * ( 2 * gap) ) + gap );
    end

    
%Calculations for Auto sorting. Also, sorts segments into the channels they
%came in as.
    ii = 1;
    for i = 1:length(breakStack)
        meanArray(ii,i-((ii-1)*divs)) = mean(cell2mat(breakStack{i}));
        stdArray(ii,i-((ii-1)*divs)) = std(cell2mat(breakStack{i}));
        sumArray(ii,i-((ii-1)*divs)) = sum(cell2mat(breakStack{i}));
        breakStackRev{ii,i-((ii-1)*divs)} = breakStack{i};
        if mod(i , divs) == 0   %breaks up Arrays into channels to simplify sorting below
            ii = ii + 1;
        end
    end

    
% Prepares temporary holding for sorting
    for i = 1:channels
        sort1{1,i} = meanArray(i,1);
        sort1{2,i} = [];
    end

%This section sorts segments into channels. Currently setup for Manual.   
totaldivs = length(breakStack);
i = 1;
while i <= totaldivs
    
    %comparing current value to each channel. For auto sorting
    for ii = 1:channels
        if i > channels
            tempDiff(ii) = abs(meanArray(i) - mean(prevBreak{ii}));
        else
            tempDiff(ii) = abs(meanArray(i) - sort1{1,ii});
        end
    end
    
    %Beginning of manual sorting
    if i > channels
        %sets location/size of figures
        figLoc = 0;
        screensize = get( groot, 'Screensize' );
        figWidth = screensize(3) / ( ceil( channels / 4 ) );
        figHeight = screensize(4) / 1.2 ; 
        figDim = [-figWidth + 1, (figHeight) - (figHeight*.9), figWidth, figHeight];
        for VV = 1:channels
            temper = cell2mat(breakStackRev{i});
            xpast = 1:1:length(prevBreak{VV});
            xpres = (1:1:length(temper))+(xpast(end));
            [~, numSections] = size(sortedCell{1,VV});
            
            %adjusts trace plot location for channel numbers greater than 4
            %4 selected because it maximizes readability and space.
            if mod(VV,4) == 1
                figLoc = figLoc + 1;
                figDim(1) = figDim(1) + figWidth ;
                figIndex = 1;
            end
            
            %Plots previous segment from each channel along with the new
            %segment in question.
            h(figLoc) = figure(555+figLoc);
            subplot(subTrace,1,figIndex)
            plot(xpast,prevBreak{VV});
            hold on
            plot(xpres,temper) 
            ylabel( num2str(numSections) )
            set(gcf, 'Position', figDim)
            figIndex = figIndex + 1;
        end
        
        % Taking user input
        title([ num2str(i) ' of ' num2str(totaldivs) ])
        prompt = 'Which channel? (Type 0 to repeat last assignment)';
        channelGuess = input(prompt);
        close(h);
        
        % Allows user to repeat previous assignment by inputting "0"
        if channelGuess == 0
            sortedCell{1,channelAssignment} = sortedCell{1,channelAssignment}(1,1:end-1);
            prevBreak = prevHold;
            i = i-1; 
            continue
        else
            % Assigns values
            channelAssignment = channelGuess;
            sortedCell{1,channelAssignment} = [sortedCell{1,channelAssignment},breakStackRev{i}];
            prevHold = prevBreak;
            prevBreak{channelAssignment} = cell2mat(breakStackRev{i});
        end
    else
        %Sorting based off closest mean value. Auto.
        channelAssignment = min(tempDiff);
        channelAssignment = find(tempDiff == channelAssignment);
        sort1{2,channelAssignment} = [sort1{2,channelAssignment},meanArray(i)];
        sort1{1,channelAssignment} = mean(sort1{2,channelAssignment});   %updates channel rep to improve sorting accuracy
        sortedCell{1,channelAssignment} = [sortedCell{1,channelAssignment},breakStackRev{i}];
        prevBreak{channelAssignment} = cell2mat(breakStackRev{i});
    end
    i = i+1;
end

%%% Preparing output data structures

% getting number of segments assigned to each channel
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

