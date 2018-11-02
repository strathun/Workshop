%% filter testing

 fs = 10e3 ;
 t = 0:( 1 / ( fs ) ):1;
 f = 1000;
 a = 4;
 y = [];
%  y(1,:)= a*sin(2*pi*f*t);
%  y(2,:)= ( a + 1) * sin(2*pi*f*t);
y = a*sin(2*pi*f*t);
z = a*2*sin(2*pi*f*t);
zz = y+z;
zz1 = zz(1:length(zz)/2);
zz2 = zz(( length(zz)/2)+ 73:end)+5;
zzz = [zz1 zz2];
y = zzz;
 %% 
 plot(y);
 ylim([-1*(a+1), a+1])
 
 %%
 
 ynew = bandpass(y.',[750 4000],fs);
 
 %%
 figure
 plot(y)
 
 figure
 plot(ynew)