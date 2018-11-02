%% This is for loading and handling Mohit's multiplexed data. 
clear

%channels of interest so far...
%ch6_2018_7_10_18_53_53_1048576_smpls_raw.mat
%ch0_1_2018_7_10_18_30_49_raw.mat
%ch5_6_2018_7_10_18_44_4_raw.mat
%2018_7_27_12_38_24_8388608_smpls_raw.mat

%ch61
%2018_8_6_13_35_59_8388608_smpls_raw.mat
%2018_8_6_13_48_36_8388608_smpls_raw.mat

% MDPI :
% 8_29 -> ch11_12
% 2018_8_29_15_35_37_8388608_smpls_raw.mat
% anaStart = 10.05;      
% anaTime = 13.95; 

% fft comparison
%2018_8_29_14_24_18_8388608_smpls_raw (2-6), (7-11)
%2018_8_29_14_24_56_8388608_smpls_raw (4-8), (9-13)
%2018_8_29_15_17_17_8388608_smpls_raw (3-7), (8-12)
%2018_8_29_15_9_43_8388608_smpls_raw (0-4), (5-9)

%
%Loads data into a usable format
V =  load('2018_8_29_15_35_37_8388608_smpls_raw.mat');  %choose channel
% V =  load('2018_8_29_14_24_18_8388608_smpls_raw.mat');  %choose channel

V1 = struct2cell(V);
V = V1{1};
V = double(V)*1.8;
V = V.';

%______User defined variables_______%
channels = 2;        % number of channels to sort
anaStart = 10.05;      %
anaTime = 13.95;       % time in seconds to analyze from recording. set to 0 to analyze full recording
% anaStart = 2;      %
% anaTime = 6;
chop = 0;            %0 for no, 1 for yes;
Fs = 600e3;          % sampling frequency
channelMax = 16;     % Based off absolute "Fs" above
threshold = -3.5;    % threshold x rms = spike detected
rejectMod = 1.7;     % reject "spikes" that are rejectMod x greater than mean shape
passBandF = 800;     % Frequency in Hz of high passband
passBandFL = 4000;   % Frequency in Hz of low passband
bpRange = [750 4000];%
bandstop = 26;       % Frequency in Hz of bandstop
order = 3;           % Order for Butterworth Filter
filterType = 4;      % 0 lowpass; 1 high pass butter; 2 bandstop; 3 bypass filter; 4 HP and LP 
avgYN = 0;           % 1 to avg data, 0 for no avg
ARP = .001;
%___________________________________%

if anaTime > 0 && anaStart > 0
    V = V( ceil( anaStart * Fs ): ceil( anaTime * Fs ) ); 
    del = mod(length(V),channels);
    V = V(1:end-del);
elseif anaTime > 0
    V = V( 1: ceil( anaTime * Fs ) ); 
    del = mod(length(V),channels);
    V = V(1:end-del);
end

% sort channels
for i = 1:channels
    Vordered(i,:) = V(i:channels:end,:);
end

%% Averaging data
% Averages data points based on the possible number of channels sampled and
% the actual number. 
if avgYN == 1
    avgPoints = channelMax / channels ; %Num points to average

    for ii = 1:channels
        Vorderedf(ii,:) = arrayfun(@(i) mean(Vordered(ii,i:i+avgPoints-1)),1:avgPoints:length(Vordered(ii,:))-avgPoints+1)'; % the averaged vector
    end
    Vordered = Vorderedf;

    Fs = Fs / channels / avgPoints ; 
else
    n = channelMax / channels ;
    Vorderedf = downsample(Vordered.',n,1);
    Vordered = Vorderedf.';
    Fs = Fs / channels / n ;
end
%% Cutting out resets and sorting
% Handles reset artefacts and sorts channels correctly.
if chop == 1
    [sortedData, breakStack, breakStackRev, eventsNew] = artefactChopper_manual(Vordered,Fs,'Manual');
    Vordered = sortedData;
elseif chop == 0
    sortedData = Vordered;
    eventsNew = [];
end

%% Filtering
% This script will filter our raw time domain data

%dataMeanSub = meanSubtraction( V );    %Removing offset, use with 3 or more channels
dataMeanSub =  Vordered;
% dataMeanSub =  V;

if filterType == 0
    %Butterworth filter _ high pass
    [ B, A ] = butter( order, passBandF / ( Fs/2 ), 'high');
    dataMeanSub = dataMeanSub.';    %filter takes column as channels, not rows
    dataHighPass = filter( B, A, dataMeanSub);
    dataHighPass = dataHighPass.';
elseif filterType == 0
    %Butterworth filter _ low pass
    [ B, A ] = butter( order, passBandFL / ( Fs/2 ));
    dataMeanSub = dataMeanSub.';    %filter takes column as channels, not rows
    dataHighPass = filter( B, A, dataMeanSub);
    dataHighPass = dataHighPass.';
elseif filterType == 2
    %Butterworth filter _ bandstop
    d = designfilt('bandstopiir','FilterOrder',6,'HalfPowerFrequency1',...
               (bandstop - 1),'HalfPowerFrequency2',(bandstop + 1), ...
               'DesignMethod','butter','SampleRate',Fs);
    dataMeanSub = dataMeanSub.';
    dataHighPass = filtfilt(d,dataMeanSub);    %bad naming convention...this isn't really high pass.
    dataHighPass = dataHighPass.';
elseif filterType == 4
    %Butterworth filter _ high pass
    [ B, A ] = butter( order, passBandF / ( Fs/2 ), 'high');
    dataMeanSub = dataMeanSub.';    %filter takes column as channels, not rows
    dataHighPass = filter( B, A, dataMeanSub);

    %Butterworth filter _ low pass
    [ B, A ] = butter( order, passBandFL / ( Fs/2 ));
    dataHighPass = filter( B, A, dataHighPass);
    dataHighPass = dataHighPass.';
elseif filterType == 5
    dataHighPass = bandpass(dataMeanSub.',bpRange,( Fs ) );
    dataHighPass = dataHighPass.';
    dataHighPass = dataHighPass(:, ( Fs*.05 ) : ( end - (Fs*.05)) );    %Cuts out filter artefacts at beginning and end of filtered data.
else
    dataHighPass = dataMeanSub ;
end

%Common Average Referencing
if channels > 3
   Vfiltered = comAvgRef( dataHighPass ); 
   %Vfiltered = ( dataHighPass ); 
else
    Vfiltered = dataHighPass;
end

%% Plots spike events

% x values for plots
% x = 1:1:length(V(1,:));
% time = (1:1:length(V(1,:)))/Fs;
x = 1:1:length(Vordered);
time = (1:1:length(Vordered))/Fs;
% timefilt = time( ( Fs*.05 ) : ( end - (Fs*.05)) );
timefilt = time;

% Detects spike events, pulls out waveforms from filtered data and plots
% along with mean waveform shape. Also plots filtered data for reference.

for i = 1:channels
    %grabbing data
    [spikesIndex, threshVal] = spike_detection(Vfiltered(i,:),threshold,1,0);
    [waveforms{i}, timeWave, spikesIndex] = waveformGrabber(Vfiltered(i,:), spikesIndex, 1.6, Fs);
    
    %spikesTemp = cell2mat(spikesIndex(i));
    spikesTemp = (spikesIndex);
    spikesTime = spikesTemp/Fs;
    spikes{i} = spikesTime;
    
    %insert chopper here
    spikesUse = spikes{1,i}.';
    [choppedData, fat] = spikeChopper(Vfiltered(i,:),spikesUse,Fs,'Threshold',1.6);
    
    Vrms(i) = rms(choppedData);
    
    figure(i)
    plot(timefilt,Vfiltered(i,:))
    hold on
    
    %adds threshold line to plot
    threshLine = ones(1,length(Vfiltered)) * threshVal;
    plot(timefilt, threshLine,'--')
    title([ 'Background Noise Vrms : ' num2str(Vrms(i)) ' uV'])
    
    %plots reset events onto graph
    for iii = 1:length(eventsNew)
        xVar = eventsNew(iii)/Fs;
        line([xVar (xVar +(xVar/1e5))],[-80 80],'Color','black')
    end
    %xlim([0 (timefilt(end))])
    xlabel('Time (s)');
    
    [events, ~] = size(waveforms{i}) ; 
    figure(i+channels)
    waveform = waveforms{i} ;
    [waveform, spikeEventsNew] = templateMatcher(waveform,rejectMod, spikesIndex, ARP, Fs); %removes "bad" spikes
    [eventsMod, ~] = size(waveform) ; 
    [threshCount(i+1), ~] = size(waveform); %+1 to add zero at beginning for below...
    spikeEventsRaster{i} = spikeEventsNew;
     
    %SNR calculation (Per RC Kelly (2007) J Neurosci 27:261)
    waveformMean = mean(waveform);
    waveformNoise = waveform - waveformMean;
    waveformNoiseSD(i) = std2(waveformNoise);
    SNRKelly(i) = ( max(waveformMean) - min(waveformMean) ) / ( 2 * waveformNoiseSD(i) );

    %SNR calculation (Per K Ludwig (2009) J Neurophys 0022-3077)
    PkPkNoiseFloor = 6 * ( std( choppedData ) );
    SNRLudwig(i) = ( max(waveformMean) - min(waveformMean) ) / PkPkNoiseFloor;
    
    for ii = 1:eventsMod
        plot(timeWave*1e3, waveform(ii,:),'Color',[.5 .5 .5], 'LineWidth', 1.2)
        hold on
    end
    threshLine = threshLine(1:(length(timeWave*1e3)));
    plot(timeWave*1e3,threshLine,'--')
    
    title([ 'SNRKelly: ' num2str(SNRKelly(i)) ', SNRLudwig: ' num2str(SNRLudwig(i)) ', SpikeCount: ' num2str(threshCount(i+1)) ])
    
    meanWave = mean(waveform) ;
    plot(timeWave*1e3, meanWave, 'LineWidth', 3.5)
    xlabel('Time (ms)');
    ylim([-80 60])
    
    figure(500)
    plot(timeWave*1e3, meanWave, 'LineWidth', 3.5)
    hold on
    xlabel('Time (ms)');
    ylim([-80 60])
%     % ISI histogram
%     isiArray(i) = diff(spikes(i))
%     figure
%     histogram(isiArray,1000)
end

% Raster Plot
figure
LineFormat.LineWidth = 0.8;
[xpoints, ypoints] = plotSpikeRaster(spikeEventsRaster,'PlotType','vertline','LineFormat',LineFormat);
hold on
    for iii = 1:length(eventsNew)
        xVar = eventsNew(iii)/Fs;
        line([xVar (xVar +(xVar/1e5))],[-80 80],'Color','blue')
    end
xlabel('Time (s)');

%% puts Raster data on filtered data graph
xpoints = (xpoints.')/Fs;
ypoints = ypoints.';
threshCount = threshCount * 3;

for i = 2: (channels+1)
    xpointsSorted = xpoints(threshCount(i-1) + 1: (threshCount(i)+ threshCount(i-1)));
    ypointsSorted = ypoints(threshCount(i-1) + 1: (threshCount(i)+ threshCount(i-1)));
    figure(i-1)
    ylim([-150 150])
    plot(xpointsSorted,((ypointsSorted)*50)-(150+(50*(i-2))),'k')
    xlim([0 4])
end
%% Plot of raw data for gut check


for i = 1:channels
    figure(300 + i)
    plot(time,Vordered(i,:))
    hold on
    xlabel('Time (s)');
end

%% Plots PSD
avgs = 64;

figure
for i = 1:channels
    [pxx1,f] = psdWalker(Vfiltered(i,:),avgs,Fs);
    loglog(f,(pxx1))
    hold on
end
ylabel('Noise Voltage ( nV / \surd Hz )','Interpreter','tex')
xlabel('Frequency (Hz)')

figure
for i = 1:channels
    [pxx1,f] = psdWalker(Vordered(i,:),avgs,Fs);
    loglog(f,(pxx1))
    hold on
end
ylabel('Noise Voltage ( nV / \surd Hz )','Interpreter','tex')
xlabel('Frequency (Hz)')


% % Plots channel raw data on the same plot. 
% figure(301)
% for i = 1:channels
%     plot(time,(1000*(i-1))+Vordered(i,:))
%     hold on
% end
% hold off

%% ISI data

for i = 1:channels
    spikeEventsDif = diff((spikeEventsRaster{i}/Fs)*(1e3));
    % [h, edges] = histcounts(spikeEventsDif,64);
    edges = 0:1:200;
    figure
    histogram(spikeEventsDif,edges)
    xlim([0 20])
    xlabel('Inter-Spike Interval (ms)')
    ylabel('Number of Counts')


    [isiArray, isiX] = isiContinuous((spikeEventsRaster{i}/Fs)*(1e3),50,4);
    figure
    plot(isiX/(1e3),isiArray)
    xlabel('Time (s)')
    ylabel('Inter-Spike Interval (ms)')
end

