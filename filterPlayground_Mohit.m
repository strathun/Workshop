% Filter testing ground. Used specifically to test filter used with Mohit's
% muxing data.

clearvars

Fs = 600e3;     %samples/s
channels = 16;  
dt = (1 / ( Fs / channels ) );
t = 0:dt:5;     % total simulation time

%Filter Specs
order = 3;
passBandF = 750;
passBandFL = 7500;

%% Generate signal
A(1,:) = 500*sin(2*pi*t*10);
A(2,:) = 200*sin(2*pi*t*100);
A(3,:) = 50*sin(2*pi*t*1000);
A(4,:) = 50*sin(2*pi*t*2000);

signal = sum(A);
%% Filter signal
%Butterworth filter _ high pass
[ B, A ] = butter( order, passBandF / ( Fs/2/16 ), 'high');
dataHighPass = filter( B, A, signal);

%Butterworth filter _ low pass
[ B, A ] = butter( order, passBandFL / ( Fs/2/16 ));
dataHighPass = filter( B, A, dataHighPass);
dataHighPass = dataHighPass.';

%%
figure
plot(t, signal)
figure
plot(t, dataHighPass)