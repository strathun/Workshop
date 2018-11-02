%Plots Impedance/Phase data similart to Gamry display.
%Right now, need to run ExtractImpedanceData2 first.

[~, ~, numTrodes] = size(Phase);

for i = 1:numTrodes
figure(i)
subplot(2,1,1)
semilogx(f(:,1,1),Phase(:,1,i)); grid on;
hold on; ylabel('Phase')
subplot(2,1,2)
loglog(f(:,1,1),Zmag(:,1,i)); grid on;
hold on; ylabel('Impedance (Magnitude)')
end