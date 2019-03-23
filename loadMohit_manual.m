%% This is for loading and handling Mohit's multiplexed data. 
clear

% channels of interest so far...
% ch6_2018_7_10_18_53_53_1048576_smpls_raw.mat
% ch0_1_2018_7_10_18_30_49_raw.mat
% ch5_6_2018_7_10_18_44_4_raw.mat
% 2018_7_27_12_38_24_8388608_smpls_raw.mat'
% 2018_7_27_12_41_11_8388608_smpls_raw.mat
% 2018_8_6_13_34_42_262144_smpls_raw.mat

%Loads data into a usable format
V =  load('2018_8_6_13_34_42_262144_smpls_raw.mat');  %choose channel
V1 = struct2cell(V);
V = V1{1};
V = double(V);
V = V.';

%______User defined variables_______%
channels = 2;        % number of channels to sort
Fs = 600000;          % sampling frequency
threshold = -4.0;    % threshold x rms = spike detected
rejectMod = 2.0;     % reject "spikes" that are rejectMod x greater than mean shape
passBandF = 750;     % Frequency in Hz of high passband
passBandFL = 6000;   % Frequency in Hz of low passband
bandstop = 26;       % Frequency in Hz of bandstop
bpRange = [750 7500];% Range for band pass filter
order = 5;           % Order for Butterworth Filter
filterType = 5;      % 0 lowpass; 1 high pass butter; 2 bandstop; 3 bypass filter; 4 HP and LP; 5 improved bandpass
sortedChannels = cell(channels,1);
%___________________________________%

% sort channels
for i = 1:channels
    Vordered(i,:) = V(i:channels:end,:);
end
% channels = 1;

% x values for plots
x = 1:1:length(Vordered);
time = (1:1:length(Vordered))/(Fs/channels);

%% Cutting into sections for manual sorting
[sortedData, breakStack, breakStackRev] = artefactChopper(Vordered,2200);

%% Filtering
% This script will filter our raw time domain data

%dataMeanSub = meanSubtraction( V );    %Removing offset, use with 3 or more channels
divisions = length(breakStack) ;

for i = 1:divisions
    
dataMeanSub =  cell2mat(breakStackRev{i}) ;

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
    dataHighPass = bandpass(dataMeanSub,bpRange,( Fs/channels ));
    dataHighPass = dataHighPass(6000:(end-6000));
else
    dataHighPass = dataMeanSub ;
end

%Common Average Referencing
if channels > 1
   Vfiltered = comAvgRef( dataHighPass ); 
   %Vfiltered = ( dataHighPass ); 
else
    Vfiltered = dataHighPass;
end

%% Plots spike events


    %grabbing data
    spikes{i} = spike_detection(Vfiltered,threshold,1);
    [waveforms{i}, timeWave] = waveformGrabber(Vfiltered, spikes{i}, 1.6, 600e3);
    
    figure(i)
    plot(Vfiltered)
    
    [events, ~] = size(waveforms{i}) ; 
    figure(i+channels)
    waveform = waveforms{i} ;
    waveform = templateMatcher(waveform,rejectMod); %removes "bad" spikes
    [eventsMod, ~] = size(waveform) ; 
    
    for ii = 1:eventsMod
        plot(timeWave*1e3, waveform(ii,:))
        hold on
    end
    
    meanWave = mean(waveform) ;
    plot(timeWave*1e3, meanWave, 'LineWidth', 3)
    xlabel('Time (ms)');
    
    prompt = 'Which channel?';
    channelGuess = input(prompt);
    sortedChannels{channelGuess} = [sortedChannels{channelGuess},Vfiltered];
    close all
end

%% Plot of raw data for gut check
% 
% figure(300)
% title('Raw Channels')
% for i = 1:channels
%     subplot(channels,1,i)
%     plot(time,Vordered(i,:))
% end
% 
% % figure(301)
% % for i = 1:channels
% %     plot(time,(1000*(i-1))+Vordered(i,:))
% %     hold on
% % end
% % hold off
