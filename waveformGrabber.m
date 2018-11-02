function [waveforms, time, spikeEventsNew] = waveformGrabber(rawdata,spikeEvents,window,Fs)
%waveformGrabber Takes raw data and spike times as input and returns waveforms of
%the spikes
%   Inputs:
%       rawdata = raw time domain data (should be filtered).
%       spikeEvents = sample number of peak waveform value.
%       window = length (in ms) of desired waveform.
%       Fs = sampling frequency.
%   Outpus:
%       waveforms : waveforms of spike events.
%%%%%
% buffer = Fs/750;
windowLength = ( ( window / 1e3 ) / ( 1 / Fs ) ) ; 
buffer = windowLength/3;
time = ( 0:1:windowLength-1 ) * (1 / Fs ) ;
[~ , events] = size( spikeEvents ) ; 
waveforms = zeros( events , windowLength );
spikeEventsNew = spikeEvents ;
% buffer = 400;
ii = 1 ; % prevents rows of zeros in waveforms.
startEvent = 1;
jj = 2;

if events > 2
    if spikeEvents( 1 ) > buffer
        beg = spikeEvents( 1 ) - buffer ;
        fin = beg + windowLength ;
        startEvent = 1;
    else
        while spikeEvents( jj ) < buffer
            jj = jj + 1;
        end
        beg = spikeEvents( jj ) - buffer ;
        fin = beg + windowLength ;
        startEvent = startEvent + 1;
        spikeEventsNew = spikeEventsNew(jj:end) ;
    end

    for i = startEvent:length( spikeEvents )
        if spikeEvents( i ) > buffer
            beg = spikeEvents( i ) - buffer ;
            fin = beg + windowLength ; 
            if fin < length(rawdata) && beg < length(rawdata)
                waveforms( ii,: ) = rawdata( beg : ( fin-1 ) ) ;
                ii = ii + 1 ; 
            end
        end

        waveforms(ii:end,:) = []; %removes excess rows of zero
    end
end

