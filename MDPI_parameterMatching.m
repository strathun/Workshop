%% Must run extractimpedancedata2 prior to this script
clear n_mins err_min alpha_min Rs_min Rp_min ii_min jj_min Hf_min

n_mins = 0;
% Bring in impedance data
ff = flipud(f(:,1,1));
Zmag = sqrt(Zreal.^2 + Zim.^2);
Zmag = flipud(Zmag(:,1,5));

Rs = 9.5e3;
Rp = 1e12;
alpha = 0.9;
err_min = inf;
for mm = 1:1 % initialize error

    % Calculate Y at 100kHz
    target = Zmag(42);
    syms Y
    eqn = target == Rs + (Rp*(1/Y)./(j*2*pi*ff(42)).^alpha)/(Rp + (1/Y)./(j*2*pi*ff(42)).^alpha);
    solY = solve(eqn,Y);
    Y = double(solY);

    % Calculate cpe values
    cpe = (1/Y)./(j*2*pi*ff).^alpha;

    % Create H(f)
    Hf = Rs + (Rp*cpe)./(Rp+cpe);
    % error
    err = sum(abs(Hf(1:42) - Zmag(1:42)));
end

n_converge = 30;
for ii = 1:1
    % increment for testing
    Rs = Rs+ii*100;
    for jj = 1:n_converge
%         Rp = 1e13;
        alpha_p = alpha + 0.02;
        alpha_m = alpha - 0.02;

        % Calculate Yp and Ym at 100kHz
        target = Zmag(42);
        syms Yp
        eqn = target == Rs + (Rp*(1/Yp)./(j*2*pi*ff(42)).^alpha)/(Rp + (1/Yp)./(j*2*pi*ff(42)).^alpha_p);
        solYp = solve(eqn,Yp);
        Yp = double(solYp);
        
        syms Ym
        eqn = target == Rs + (Rp*(1/Ym)./(j*2*pi*ff(42)).^alpha)/(Rp + (1/Ym)./(j*2*pi*ff(42)).^alpha_m);
        solYm = solve(eqn,Ym);
        Ym = double(solYm);

        % Calculate cpe values
        cpe_p = (1/Yp)./(j*2*pi*ff).^alpha_p;
        cpe_m = (1/Ym)./(j*2*pi*ff).^alpha_m;
        
        % Create H(f)
        Hf_p = Rs + (Rp*cpe_p)./(Rp+cpe_p);
        Hf_m = Rs + (Rp*cpe_m)./(Rp+cpe_m);
        
        % error
        err_p = sum(abs(Hf_p(1:42) - Zmag(1:42)));
        err_m = sum(abs(Hf_m(1:42) - Zmag(1:42)));
        
        if err_p > err_m
            alpha = alpha_p;
            err = err_p;
        elseif err_m > err_p
            alpha = alpha_m;
            err = err_m;
        end
        
        if err < err_min
            n_mins = n_mins+1;
            err_min(n_mins) = err;
            alpha_min(n_mins) = alpha;
            Rs_min(n_mins) = Rs;
            Rp_min(n_mins) = Rp;
            ii_min(n_mins) = ii;
            jj_min(n_mins) = jj;
            Hf_min(:,n_mins) = Hf;
            
            if jj > n_converge - 10
                figure()
                loglog(ff,Hf_p)
                hold on;
                loglog(ff,Zmag)
            end
        end
        

    end
    ii
end


loglog(ff,Hf_p)
hold on;
loglog(ff,Zmag)


