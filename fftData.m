% fft for time domain

%% raw data
%load data files
% load('')
% load('')
% 
% %channel selection
% V = rawdata( 2,:);
% V = double( V ) / ( ( 1e6 ) * 4 );

% V = Vfiltered( 3, 1800100*2:1800100*3);
% V = double( V ) / ( ( 1e6 ) * 4 );

%for comAvgRef filtered data
V = Vfiltered( 4,: );
% % V = quiet( 1,: );
% V = V/( 4 * ( 1e6 ) );

Vwelch = V;
V = V.^2;
% figure
% plot(time,V); 

%% fft

Fs = 30000;                 % Sampling frequency                    
T = 1/Fs;                   % Sampling period       
L = length(V);              % Length of signal
windowTime = 10; %seconds
windowLength = windowTime * Fs;
windowNum = floor( L / ( windowLength ) );

freq = Fs * ( 0:( L/2 ) ) / L;
df = freq( 2 );
freq = windowNum*freq(1:windowLength/2);

for win = 1:windowNum
    ffV(win,:) = fft(V( (((win-1)*windowLength)+1):(win*windowLength)));
    P2 = abs(ffV(win,:)/windowLength);
    P1 = P2(1:windowLength/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    ffVstack(win,:) = P1;
end

ffVavg = mean(ffVstack);
ffVavg = sqrt(ffVavg(2:end));

% f = (0:L-1)/(L*T);
% df = f(2);
% fftV = fft(V)/L/sqrt(df);
% PSDV = (fftV.*conj(fftV));

figure
subplot(2,2,1)
semilogx(freq,ffVavg) 
xlim([100 15e3])

%% fft test using autocorr method

Fs = 30000;
t = 0:1/Fs:1-1/Fs;

Rxx = xcorr(V,'biased');
Rxxdft = sqrt(abs(fftshift(fft(Rxx))));
freq = -Fs/2:Fs/length(Rxx):Fs/2-(Fs/length(Rxx));

subplot(2,2,2)
semilogx(freq,Rxxdft);
xlim([100 15e3])

%% power using Welch method
% use this one!!!
load('groundPower');
%for i = 1:16
% x = V(4,(1800e3)*1:(1800e3)*2);
% na = 64; %number of averages
% nx = max(size(x));
% w = floor(nx/na);
% [pxx, f] = pwelch(x,w,[],[],30000);

%[pxx, f] = pwelch(chopedData( 1,:),8193,[],[],30000);
% [pxx, f] = pwelch(V( 4,(1800e3)*2:end),8193,[],[],30000);
[pxx, f] = pwelch(V( 6,(900e3)*5:end),8193,[],[],30000);

pxx = pxx - groundPower;
pxx1 = (sqrt(pxx))*1e9; %took out the /2 for now... Not really sure what the correct method is.

%  noiseArray(i,:) = pxx1 ;
%  powerArray(i,:) = pxx ; 
%  avgNoise = mean(noiseArray);
%  verifArray(1,1) = trapz(pxx(176:21847)*f(2));
% verifArray(i,2) = var(V(i,:));
% verifArray(i,3) = (abs(verifArray(i,1)-verifArray(i,2))/verifArray(i,2))*100;

%figure(i)
figure(2)
%subplot( 2, 1, 1 )
loglog(f,(pxx1))
xlim([0 15e3])
ylim([1 100e3])
ylabel('Noise Voltage ( nV / \surd Hz )','Interpreter','tex')
xlabel('Frequency (Hz)')
legend('Still', 'Breathing')
hold on
%end

%%
%Use the next 2 lines for testing the data processing. Comment when running
%full script.
s.Rate = 30000;
s.DurationInSeconds = 60*3;

%NFFT = 2^(nextpow2(s.Rate));
NFFT = s.Rate;

%Generate power spectrum for signal

Y = zeros(s.DurationInSeconds, NFFT/2);
Y1 = zeros(1, NFFT);

freqScale = (s.Rate)*(0:(NFFT/2))/NFFT;

for i = 1:s.DurationInSeconds

Y1 = fft((V(((i-1)*s.Rate)+1:(i*s.Rate))),NFFT); %1/4
P2 = abs(Y1);
% P2 = abs(Y1);
P1 = P2(1:(NFFT/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
P1 = P1(2:end);
Y(i,:) = P1;

end
Y = mean(Y);
%taking the average of each second 
%Y = mean(Y)/sqrt(499e3)/380;

% Do power/bin to V conversion
% Y = Y/sqrt(499e3);   %1/4
freqScale = freqScale(2:end);

%just for test. Remove during full script run
subplot(2,2,4)
semilogx(freqScale, abs(Y))
ylabel('V_{rms}')
xlim([100 15e3])

%%
[Pxx,F] = periodogram(sqrt(V),[],length(V),Fs);
figure;
semilogx(F,Pxx)
xlim([100 15e3])

%%
psd = abs(fft(sqrt(V))).^2/length(V);

figure
freq = -Fs/2:Fs/length(psd):Fs/2-(Fs/length(psd));
semilogx(freq(1:length(freq)-1),psd(1:length(freq)-1))


