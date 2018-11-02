function [sortedArray, spikeEventsNew] = templateMatcher(array, devs, spikeEvents, ARP, Fs)
%UNTITLED2 Summary of this function goes here
%   Inputs : 
%       array = channels should be arranged as rows
%       devs = # of std deviations away to be rejected
%       spikeEvents = 
%       ARP = absolute refractory period. Removes any waveforms occuring
%           less than this value (seconds).
%%%

[traces, points] = size(array);

    if traces > 1
        meanArray = mean(array) ; 
    else
        meanArray = array ;
    end
    
devArray = ( std(array) ) * devs;
highArray = meanArray + devArray ; 
lowArray = meanArray - devArray ;

sortedArray = array; 

for i = 1:traces
    ii = 1 ; 
    while ii < points + 1
        if array(i,ii) > highArray(ii) || array(i,ii) < lowArray(ii)
            sortedArray(i,:) = 0 ;
            spikeEvents(i) = 0 ;
            ii = points + 2;
        end
        ii = ii + 1;
    end
end

%removes waveforms that occur before ARP
spikeTimes = spikeEvents / Fs ;
isiArray = diff(spikeTimes);
throwAway = find(isiArray <= ARP);
if size(throwAway) > 0
    for i = 1:length(throwAway)
        eventIndex = throwAway(i) + 1 ;
        sortedArray(eventIndex,:) = 0 ;
        spikeEvents(eventIndex) = 0 ;
    end
end

spikeEventsNew = spikeEvents ;
spikeEventsNew(spikeEventsNew==0) = [];
sortedArray = sortedArray(any(sortedArray,2),:);
end

