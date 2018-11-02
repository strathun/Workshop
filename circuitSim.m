f = logspace(0,5,1000);
j = sqrt(-1);


R1 = 4.6e3 ;
R2 = 3.5e6;
C = .03e-6;

for i = 1:length(f)
    s = j*2*pi*f(i);
    A = ( R1*R2*C*s ) + R1 + R2;
    B = ( R2*C*s ) + 1;
    imp(i) = A/B;
end
loglog(f,imp);
% Zmag = 