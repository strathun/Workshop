%

figure(6)
fig = gcf;
fig = get( fig, 'Number' );
for ii = 1:4
    startTime = (ii-1)*7;
    stopTime = ii*7;
%     fname1 = sprintf('Ch.%02d_time_%02d-%02d.fig',chSelect ,startTime, stopTime);
    fname2 = sprintf('Ch.%02d_time_%02d-%02d.png',chSelect ,startTime, stopTime);
    xlim([startTime stopTime])
%     saveas(figure(fig), fname1)
    saveas(figure(fig), fname2)
end
