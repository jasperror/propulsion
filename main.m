% ENGR477 Propulsion
% Project 1: F135 Engine
% Jasper Palmer & Jackson Gilbert

clear; close all; clc;

%% Constants

% Design-point Operating Conditions
% dmdt_a = 150; % Total air mass flow rate, kg/s
% pr = 28; % Overall pressure ratio, []
% T_5 = 2000; % TIT, K
% BPR = 0.57; % Bypass ratio, []
FAN_pr = 1.75; % Fan pressure ratio, []
% LPC_pr = 1.25; % LPC pressure ratio, []
% HPC_pr = 12.8; % HPC pressure ratio, []
LHV = 43150; % Fuel heating value, kJ/kg
% T_a = 288.15; % Design point ambient temperature, K
% p_a = 101.325e3; % Ambient pressure, Pa

% Flow & Geometric Assumptions
% M_1 = 0.5; % Fan inlet mach number, []
% Vbar_45 = 150; % Combustor average axial velocity, m/s
d_9 = 0.78; % Nozzle throat diameter, m
d_10 = 0.78; % Nozzle exit diameter, m
% M_7 = 0.5; % Turbine exit mach number, []

% Component Efficiencies
eta.FAN = 0.89; % Fan efficiency, []
eta.LPC = 0.88; % LPC efficiency, []
eta.HPC = 0.86; % HPC efficiency, []
eta.BRN = 0.99; % Burner efficiency, []
eta.HPT = 0.89; % HPT efficiency, []
eta.LPT = 0.91; % LPT efficiency, []
eta.NOZ = 0.98; % Nozzle efficiency, []
eta.NOZWET = 0.97; % Nozzle efficiency with afterburner, []
eta.SFT = 0.99; % Shaft efficiency, []

% Stagnation Pressure Ratios
spr.INT = 0.99; % Intake SPR, []
spr.LPC = 0.99; % LPC duct SPR, []
spr.BRN = 0.94; % Burner SPR, []
spr.MXR = 0.97; % Mixer SPR, []
spr.BPD = 0.96; % Bypass duct SPR, []
spr.HPC = 0.99; % HPC duct SPR, []
spr.ABR = 0.98; % Afterburner duct SPR, []
spr.ABRON = 0.95; % Afterburner operational SPR, []
spr.NOZ = 0.98; % Nozzle SPR, []

%% Engine Components

function out = diffuser(in, spr, Ar, R) % Assuming adiabatic
    out.p0 = in.p0 * spr;
    out.dmdt = in.dmdt;
    out.T0 = in.T0;

    i = 0; err = 1; tmp = in.T; tmpp = in.p;
    while err > 0.001 && i < 1e4
        i = i+1;
        cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',tmp,'P',tmpp,'Air');
        cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',tmp,'P',tmpp,'Air');
        out.gamma = cp/cv;
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',tmp,'P',tmpp,'Air');
        out.V = in.rho/out.rho*Ar*in.V;
        out.M = out.V/sqrt(out.gamma*R*tmp);
        out.p = out.p0*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(1-out.gamma));
        out.T = in.T0/(1+(out.gamma-1)/2*out.M^2);
        err = sqrt(mean( ((out.T-tmp)/out.T)^2+((out.p-tmpp)/out.p)^2 ));
        tmp = out.T;
        tmpp = out.p;
    end 
    out.h = py.CoolProp.CoolProp.PropsSI('H','T',out.T,'P',out.p,'Air');
end

function out = comp(in, pr, eta, Ar, R)
    out.p = in.p*pr;
    out.dmdt = in.dmdt;
    h1 = py.CoolProp.CoolProp.PropsSI('H','T',in.T,'P',in.p,'Air');
    h01 = h1 + in.V^2/2;
    s1 = py.CoolProp.CoolProp.PropsSI('S','T',in.T,'P',in.p,'Air');
    
    tmp = in.T; i = 0; err = 1;
    while err > 0.001 && i < 1e4
        i = i + 1;
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',tmp,'P',out.p,'Air');
        out.V = in.rho/out.rho*Ar*in.V;
        h2s = py.CoolProp.CoolProp.PropsSI('H','P',out.p,'S',s1,'Air');
        h02s = h2s + out.V^2/2;
        h02 = (h02s-h01)/eta+h01;
        h2 = h02 - out.V^2/2;
        out.T = py.CoolProp.CoolProp.PropsSI('T','H',h2,'P',out.p,'Air');
        err = abs((out.T-tmp)/out.T);
        tmp = out.T;
    end
    
    cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',out.T,'P',out.p,'Air');
    cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',out.T,'P',out.p,'Air');
    out.gamma = cp/cv;
    out.h = h2;
    out.h0 = h02;
    out.M = out.V/sqrt(out.gamma*R*out.T);
    if out.M > 1
        warning('Mach>1 in compressor stage!')
    end
    out.p0 = out.p*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.T0 = out.T*(1+(out.gamma-1)/2*out.M^2);
    out.w = h2-h1+0.5*(out.V^2-in.V^2);
end

function out = turb(in, w_req, eta, Ar, R)
    out.dmdt = in.dmdt;
    h1 = py.CoolProp.CoolProp.PropsSI('H','T',in.T,'P',in.p,'Air');
    h01 = h1 + in.V^2/2;
    s1 = py.CoolProp.CoolProp.PropsSI('S','T',in.T,'P',in.p,'Air');
    h02 = h01-w_req;
    h02s = h01-(h01-h02)/eta;
    
    tmp = in.T; tmpp = in.p; i = 0; err = 1;
    while err > 0.001 && i < 1e4
        i = i + 1;
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',tmp,'P',tmpp,'Air');
        out.V = in.rho/out.rho*Ar*in.V;
        h2s = h02s-out.V^2/2;
        h2 = h02-out.V^2/2;
        out.p = py.CoolProp.CoolProp.PropsSI('P','H',h2s,'S',s1,'Air');
        out.T = py.CoolProp.CoolProp.PropsSI('T','H',h2,'P',tmpp,'Air');
        err = sqrt(mean( ((out.p-tmpp)/out.p)^2+((out.T-tmp)/out.T)^2 ));
        tmp = out.T;
        tmpp = out.p;
    end

    cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',out.T,'P',out.p,'Air');
    cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',out.T,'P',out.p,'Air');
    out.gamma = cp/cv;
    out.M = out.V/sqrt(out.gamma*R*out.T);
    out.p0 = out.p*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.T0 = out.T*(1+(out.gamma-1)/2*out.M^2);
    out.h = h2;
    out.h0 = h02;
end

function out = duct(in, spr, R)
    out.dmdt = in.dmdt;
    out.p0 = in.p0 * spr;
    out.T0 = in.T0;
    tmp = in.T; tmpp = in.p; i = 0; err = 1;
    while err > 0.001 && i < 1e4
        i = i + 1;
        cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',tmp,'P',tmpp,'Air');
        cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',tmp,'P',tmpp,'Air');
        out.gamma = cp/cv;
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',tmp,'P',tmpp,'Air');
        out.V = in.rho/out.rho*in.V; % Assuming no change in cross-section area
        out.M = in.V/sqrt(out.gamma*R*tmp);
        out.p = out.p0*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(1-out.gamma));
        out.T = out.T0/(1+(out.gamma-1)/2*out.M^2);
        err = sqrt(mean( ((out.T-tmp)/out.T)^2+((out.p-tmpp)/out.p)^2 ));
        out.h = py.CoolProp.CoolProp.PropsSI('H','T',out.T,'P',out.p,'Air');
        out.h0 = out.h + out.V^2/2;
        tmp = out.T;
        tmpp = out.p;
    end
end

function out = CEARUN(p, T, CEA, Tv)
    % Import CEA data
    sz = size(CEA, [1, 3]);
    Tad = reshape(CEA(:,4,:),sz);
    c = reshape(CEA(:,7,:),sz);
    h_BNR = reshape(CEA(:,5,:),sz);
    gamma = reshape(CEA(:,6,:),sz);
    CO = reshape(CEA(:,8,:),sz);
    CO2 = reshape(CEA(:,9,:),sz);
    NO = reshape(CEA(:,10,:),sz);
    OF = reshape(CEA(1,2,:),[1, size(CEA,3)]);
    phi = reshape(CEA(1,1,:),[1, size(CEA,3)]);
    pv = CEA(:,3,1)';
    h_air = py.CoolProp.CoolProp.PropsSI('H','T',T,'P',p,'Air')/1000-py.CoolProp.CoolProp.PropsSI('H','T',298.15,'P',101325,'Air')/1000; % Air absolute enthalpy [kJ/kg]
    h_fuel = -284117/151.9; % Jet-A absolute enthalpy [kJ/kg]

    % Lookup Values
    if T < 200 % If outside of bounds then use min/max. There is probably a better way to do this
        out.OF = CEA(1,2,1);
        out.phi = CEA(1,1,1);
        h_reac = (h_air*OF + h_fuel)/(1+OF);
        if p*1e-5 < 0.1
            out.T = CEA(1,4,1);
            out.c = CEA(1, 7, 1);
            out.h = CEA(1, 5, 1);
            out.Q_R = out.h-h_reac;
            out.gamma = CEA(1,6,1);
            out.CO2e = 10*CEA(1,8,1)+CEA(1,9,1);
            out.NO = CEA(1,10,1);
        elseif p*1e-5 > 50
            out.T = CEA(end,4,1);
            out.c = CEA(end, 7, 1);
            out.h = CEA(end, 5, 1);
            out.Q_R = out.h-h_reac;
            out.gamma = CEA(end,6,1);
            out.CO2e = 10*CEA(end,8,1)+CEA(end,9,1);
            out.NO = CEA(end,10,1);
        else
            out.T = interp1(pv, CEA(:,4,1), p*1e-5);
            out.c = interp1(pv, CEA(:,7,1), p*1e-5);
            out.h = interp1(pv, CEA(:,5,1), p*1e-5);
            out.Q_R = out.h-h_reac;
            out.gamma = interp1(pv, CEA(:,6,1), p*1e-5);
            out.CO2e = 10*interp1(pv, CEA(:,8,1), p*1e-5)+interp1(pv, CEA(:,9,1), p*1e-5);
            out.NO = interp1(pv, CEA(:,10,1), p*1e-5);
        end
    elseif T > 2000
        out.OF = CEA(1,2,end);
        out.phi = CEA(1,1,end);
        h_reac = (h_air*OF + h_fuel)/(1+OF);
        if p*1e-5 < 0.1
            out.T = CEA(1,4,end);
            out.c = CEA(1, 7, end);
            out.h = CEA(1, 5, end);
            out.Q_R = out.h-h_reac;
            out.gamma = CEA(1,6,end);
            out.CO2e = 10*CEA(1,8,end)+CEA(1,9,end);
            out.NO = CEA(1,10,end);
        elseif p*1e-5 > 50
            out.T = CEA(end,4,end);
            out.c = CEA(end, 7, end);
            out.h = CEA(end, 5, end);
            out.Q_R = out.h-h_reac;
            out.gamma = CEA(end,6,end);
            out.CO2e = 10*CEA(end,8,end)+CEA(end,9,end);
            out.NO = CEA(end,10,end);
        else
            out.T = interp1(pv, CEA(:, 4, end), p*1e-5);
            out.c = interp1(pv, CEA(:, 7, end), p*1e-5);
            out.h = interp1(pv, CEA(:, 5, end), p*1e-5);
            out.Q_R = out.h-h_reac;
            out.gamma = interp1(pv, CEA(:,6,end), p*1e-5);
            out.CO2e = 10*interp1(pv, CEA(:,8,end), p*1e-5)+interp1(pv, CEA(:,9,end), p*1e-5);
            out.NO = interp1(pv, CEA(:,10,end), p*1e-5);
        end
    else
        out.OF = interp1(Tv, OF, T);
        out.phi = interp1(Tv, phi, T);
        h_reac = (h_air*OF + h_fuel)/(1+OF);
        if p*1e-5 < 0.1
            size(Tv)
            size(CEA(1,4,:))
            out.T = interp1(Tv, reshape(CEA(1,4,:), [1, size(CEA,3)]), T);
            out.c = interp1(Tv, reshape(CEA(1,7,:), [1, size(CEA,3)]), T);
            out.h = interp1(Tv, reshape(CEA(1,5,:), [1, size(CEA,3)]), T);
            out.Q_R = out.h-h_reac;
            out.gamma = interp1(Tv, reshape(CEA(1,6,:), [1, size(CEA,3)]), T);
            out.CO2e = 10*interp1(Tv, reshape(CEA(1,8,:), [1, size(CEA,3)]), T)+interp1(Tv, reshape(CEA(1,9,:), [1, size(CEA,3)]), T);
            out.NO = interp1(Tv, reshape(CEA(1,10,:), [1, size(CEA,3)]), T);
        elseif p*1e-5 > 50
            out.T = interp1(Tv, reshape(CEA(end,4,:), [1, size(CEA,3)]), T);
            out.c = interp1(Tv, reshape(CEA(end,7,:), [1, size(CEA,3)]), T);
            out.h = interp1(Tv, reshape(CEA(end,5,:), [1, size(CEA,3)]), T);
            out.Q_R = out.h-h_reac;
            out.gamma = interp1(Tv, reshape(CEA(end,6,:), [1, size(CEA,3)]), T);
            out.CO2e = 10*interp1(Tv, reshape(CEA(end,8,:), [1, size(CEA,3)]), T)+interp1(Tv, reshape(CEA(end,9,:), [1, size(CEA,3)]), T);
            out.NO = interp1(Tv, reshape(CEA(end,10,:), [1, size(CEA,3)]), T);
        else % If parameters are OK then interpolate values from CEA
            [X, Y] = meshgrid(Tv, pv);
            out.T = interp2(X, Y, Tad, T, p*1e-5);
            out.c = interp2(X, Y, c, T, p*1e-5);
            out.h = interp2(X, Y, h_BNR, T, p*1e-5);
            out.Q_R = out.h-h_reac;
            out.gamma = interp2(X, Y, gamma, T, p*1e-5);
            out.CO2e = 10*interp2(X, Y, CO, T, p*1e-5) + interp2(X, Y, CO2, T, p*1e-5);
            out.NO = interp2(X, Y, NO, T, p*1e-5);
        end
    end
    out.Q_R = -(py.CoolProp.CoolProp.PropsSI('H','T',out.T,'P',p,'Air')/1000-py.CoolProp.CoolProp.PropsSI('H','T',T,'P',p,'Air')/1000); % super scuffed not using the CEA (impove me!)
end

function out = combustor(in, spr, LHV, eta, R, Tad)
    i = 0; err = 1; tmp = in.T0; tmpp = in.p0;
    while err > 0.001 && i < 1e4
        % disp(['i=',num2str(i),', T=',num2str(tmp)])
        i = i+1;
        cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',tmp,'P',tmpp,'Air');
        cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',tmp,'P',tmpp,'Air');
        gamma = cp/cv;
        M = in.V/sqrt(gamma*R*tmp);
        T = in.T0/(1+(gamma-1)/2*M^2);
        p = in.p0*(1+(gamma-1)/2*M^2)^(-gamma/(gamma-1));
        err = sqrt(mean( ((T-tmp)/T)^2+((p-tmpp)/p)^2 ));
        tmp = T;
        tmpp = p;
    end
    
    CEA = load('CEA.mat');
    switch Tad
        case 1750
            CEA = CEA.COMB175;
            Tv = [200, 800, 1000, 1500];
        case 2000
            CEA = CEA.COMB;
            Tv = [200, 500, 700, 800, 900, 1000, 2000];
        case 2250
            CEA = CEA.COMB225;
            Tv = [200, 800, 1000, 1500];
    end
    out = CEARUN(p, T, CEA, Tv);
    out.p = p;
    out.h = py.CoolProp.CoolProp.PropsSI('H','T',out.T,'P',out.p,'Air');
    out.p0 = spr*in.p0;
    out.rho = py.CoolProp.CoolProp.PropsSI('D','T',out.T,'P',out.p,'Air');

    % Calculate resulting parameters
    out.V = in.V;
    out.M = out.V/out.c;
    t0_star = in.T0*(1+out.gamma*out.M^2)^2/(2*(out.gamma+1)*out.M^2*(1+(out.gamma-1)/2*out.M^2));
    out.T0 = t0_star*(in.T0/t0_star + out.Q_R/(out.gamma*R/(out.gamma-1)*t0_star));
    out.dmdt_f = -out.Q_R*in.dmdt/LHV/eta;
    out.dmdt = out.dmdt_f+in.dmdt;
    out.h0 = out.h + out.V^2/2;
end

function out = mix(core, bypass, spr, R)
    hb = py.CoolProp.CoolProp.PropsSI('H','T',bypass.T,'P',bypass.p,'Air');
    hc = py.CoolProp.CoolProp.PropsSI('H','T',core.T,'P',core.p,'Air');
    h0b = hb+bypass.V^2/2;
    h0c = hc+core.V^2/2;
    h02 = (bypass.dmdt*h0b+core.dmdt*h0c)/(bypass.dmdt+core.dmdt);

    d = 1;
    r = d/2;
    AH = pi*(0.85*r)^2; % Assuming core outlet is 85% of radius of engine
    A2 = pi*r^2;
    AC = A2-AH;
    out.dmdt = core.dmdt+bypass.dmdt;
    
    tmp=(core.rho+bypass.rho)/2; i=0; err=1;
    while err > 0.001 && i < 1e4
        i = i + 1;
        out.V = (bypass.rho*AC*bypass.V+core.rho*AH*core.V)/(tmp*A2);
        h2 = h02-out.V^2/2;
        out.p = (out.dmdt*out.V-bypass.dmdt*bypass.V-core.dmdt*core.V+bypass.p*AC+core.p*AH)/A2;
        out.rho = py.CoolProp.CoolProp.PropsSI('D','P',out.p,'H',h2,'Air');
        err = abs((out.rho-tmp)/out.rho);
        tmp = out.rho;
    end

    cp = py.CoolProp.CoolProp.PropsSI('CPMASS','H',h2,'P',out.p,'Air');
    cv = py.CoolProp.CoolProp.PropsSI('CVMASS','H',h2,'P',out.p,'Air');
    out.gamma = cp/cv;
    out.T = py.CoolProp.CoolProp.PropsSI('T','H',h2,'P',out.p,'Air');
    out.M = out.V/sqrt(out.gamma*R*out.T);
    out.p0 = spr*out.p*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.T0 = out.T*(1+(out.gamma-1)/2*out.M^2);
    out.h = h2;
end

function out = afterburner_inop(in, spr, R)
    out = duct(in, spr, R);
    out.dmdt_f = 0;
    out.Q_R = 0;
end

function out = afterburner_op(in, spr, LHV)
    CEA = load('CEA.mat');
    out = CEARUN(in.p, in.T, CEA.ABR, [200, 500, 700, 900, 2000]);
    out.p0 = in.p0*spr;
    out.dmdt_f = -out.Q_R*in.dmdt/LHV;
    out.dmdt = in.dmdt+out.dmdt_f;
    i = 0; err = 1; tmp = out.p0;
    while err > 0.001 && i < 1e4
        i = i + 1;
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',out.T,'P',tmp,'Air');
        out.V = in.rho/out.rho*in.V;
        out.M = out.V/out.c;
        out.T0 = out.T*(1+(out.gamma-1)/2*out.M^2);
        out.p = out.p0*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(1-out.gamma));
        err = abs((out.p-tmp)/out.p);
        tmp = out.p;
    end
end

function out = nozzle(in, eta, A, R, p_amb)
    
    out = in;
    mdot = in.dmdt;
    h1 = py.CoolProp.CoolProp.PropsSI('H','T',in.T,'P',in.p,'Air');
    h01 = h1 + in.V^2/2;
    p01 = in.p0;
    T01 = in.T0;
    gamma = in.gamma;
    pr_crit = (2/(gamma+1))^(gamma/(gamma-1));
    
    if p_amb/p01 <= pr_crit
        
        %Choked Flow
    
        out.M = 1;
        out.T = T01/(1+(gamma-1)/2);
        out.p = p01*(2/(gamma+1))^(gamma/(gamma-1));
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',out.T,'P',out.p,'Air');
        a = sqrt(gamma*R*out.T);
        out.V = a;
        out.h = py.CoolProp.CoolProp.PropsSI('H','T',in.T,'P',in.p,'Air');
        out.h0 = out.h + a^2/2;
    
    else

    %Unchhoked Flow
    
        out.p = p_amb;
        Tguess = in.T*0.8;
        err = 1;
        i = 0;
    
        while err > 1e-5 && i < 1000
    
            i = i + 1;
            rho = py.CoolProp.CoolProp.PropsSI('D','T',Tguess,'P',out.p,'Air');
            V = mdot/(rho*A);
            h2 = h01 - V^2/(2*eta);
            Tnew = py.CoolProp.CoolProp.PropsSI('T','H',h2,'P',out.p,'Air');
            err = abs((Tnew-Tguess)/Tnew);
            Tguess = Tnew;
    
        end
    
        out.T = Tnew;
        out.V = V;
        out.rho = rho;
        cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',out.T,'P',out.p,'Air');
        cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',out.T,'P',out.p,'Air');
        out.gamma = cp/cv;
        out.M = out.V/sqrt(out.gamma*R*out.T);
        out.h = h2;

    end

    out.p0 = out.p*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.T0 = out.T*(1+(out.gamma-1)/2*out.M^2);

end

%% Engine Structure

%{
Element Numbering
a: Atmospheric (inlet condition)
1: Diffuser outlet
2: Fan outlet
3: LPC outlet
4: HPC outlet
5: Combustor outlet
6: HPT outlet
7: Bypass outlet
7': LPT outlet
8: Mixer outlet
9: Afterburner outlet
10: Nozzle outlet
%}

%{
Shafts
Shaft 1: LPC, LPT, Fan
Shaft 2: HPC, HPT
%}

%{
Afterburning
Total airflow increases to 165 kg/s
ABR reheats flow to 2450K
Nozzle throat diameter increased to 0.92m and exit diameter to 1.15m
Afterburner total pressure ratio is reduced to 0.95
Nozzle efficiency reduced to 0.97
All other design point parameters, efficiencies, and pressure-loss
assumptions are identical to those listed for dry operations.
%}

%{
Fuel type: Jet-A/JP-8
%}

%{
% Ambient Conditions
ambient.T = 15+273.15;
ambient.p = 101.325e3;
ambient.dmdt = 150;
ambient.rho = py.CoolProp.CoolProp.PropsSI('D','T',ambient.T,'P',ambient.p,'Air');
R = 287;
ambient.M = 0.5;
cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',ambient.T,'P',ambient.p,'Air');
cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',ambient.T,'P',ambient.p,'Air');
ambient.h = py.CoolProp.CoolProp.PropsSI('H','T',ambient.T,'P',ambient.p,'Air');
ambient.gamma = cp/cv;
ambient.V = ambient.M*sqrt(ambient.gamma*R*ambient.T);
ambient.h0 = ambient.h + ambient.V^2/2;
dmdt_aH = ambient.dmdt/(BPR+1);
dmdt_aC = BPR*ambient.dmdt/(BPR+1);
u = ambient.M*sqrt(ambient.gamma*R*ambient.T);
ambient.p0 = ambient.p*(1+(ambient.gamma-1)/2*ambient.M^2)^(ambient.gamma/(ambient.gamma-1));
ambient.T0 = ambient.T*(1+(ambient.gamma-1)/2*ambient.M^2);

% Engine Calculations
engine.diffuser = diffuser(ambient, spr.INT, 0.95, R);
engine.fan = comp(engine.diffuser, FAN_pr, eta.FAN, 1.02, R);

core = engine.fan;
core.dmdt = dmdt_aH;
core = duct(core, spr.LPC, R);
bypass = engine.fan;
bypass.dmdt = dmdt_aC;
bypass = duct(bypass, spr.BPD, R);


engine.lpc = comp(core, LPC_pr, eta.LPC, 3, R);
engine.lpc_ducted = duct(engine.lpc, spr.HPC, R);
engine.hpc = comp(engine.lpc_ducted, HPC_pr, eta.HPC, 3, R); % Adjusted to get 150m/s in combustor in static conditions
engine.combustor = combustor(engine.hpc, spr.BRN, LHV, eta.BRN, R, 2000);
engine.hpt = turb(engine.combustor, engine.hpc.w/eta.SFT, eta.HPT, 0.35, R);
engine.lpt = turb(engine.hpt, (engine.lpc.w*dmdt_aH+engine.fan.w*ambient.dmdt)/dmdt_aH/eta.SFT, eta.LPT, 0.7, R);
engine.mixer = mix(engine.lpt, bypass, spr.MXR, R);

engine.afterburner = afterburner_inop(engine.mixer, spr.ABR, R);
engine.afterburner.Q_R = 0;
engine.afterburner_op = afterburner_op(engine.mixer, spr.ABRON, LHV);
engine.nozzle = nozzle(engine.afterburner, eta.NOZ, pi*(0.78/2)^2, R, ambient.p);
engine.nozzle_op = nozzle(engine.afterburner_op, 0.97, pi*(0.92/2)^2, R, ambient.p);

% Output Parameters

thrust = (ambient.dmdt+engine.combustor.dmdt_f+engine.afterburner.dmdt_f)*engine.nozzle.V-ambient.dmdt*ambient.V+pi*d_10^2/4*(engine.nozzle.p-ambient.p); % Thrust, N
ST = thrust/ambient.dmdt; % Specific thrust, Ns/kg
%TSFC = (engine.combustor.dmdt_f+engine.afterburner.dmdt_f)/thrust; % Thrust-specific fuel consumption, kg/Ns
mdot_total = ambient.dmdt + engine.combustor.dmdt_f + engine.afterburner.dmdt_f;
mdot_f_total = engine.combustor.dmdt_f + engine.afterburner.dmdt_f;

f = (engine.combustor.dmdt_f+engine.afterburner.dmdt_f)/ambient.dmdt;
eta_th = (mdot_total*(engine.nozzle.V^2/2 - ambient.V^2/2)/(mdot_f_total*LHV*1000));
eta_p = thrust*ambient.V/ambient.dmdt/((1+f)*(engine.nozzle.V^2/2)-ambient.V^2/2); % Propulsion efficiency, []
eta_0 = eta_th*eta_p; % Overall efficiency, []

g=9.81;
LD = 10; % Estimate of cruise lift to drag ratio, []
s = -(eta_0*LD*engine.combustor.Q_R/g*log(1.66)); % Aircraft range (assuming fuel is 40% of aircraft weight), m

thrust = [
    (ambient.dmdt+engine.combustor.dmdt_f)*engine.nozzle.V-ambient.dmdt*ambient.V+pi*d_10^2/4*(engine.nozzle.p-ambient.p);
    (ambient.dmdt+engine.combustor.dmdt_f+engine.afterburner_op.dmdt_f)*engine.nozzle.V-ambient.dmdt*ambient.V+pi*1.15^2/4*(engine.nozzle.p-ambient.p);
]; % Thrust, N
ST = thrust./ambient.dmdt; % Specific thrust, Ns/kg
TSFC = [engine.combustor.dmdt_f; (engine.combustor.dmdt_f+engine.afterburner_op.dmdt_f)]./thrust; % Thrust-specific fuel consumption, kg/Ns

DKE = [
    (ambient.dmdt+engine.combustor.dmdt_f)*engine.nozzle.V^2/2-ambient.dmdt*ambient.V^2/2;
    (ambient.dmdt+engine.afterburner_op.dmdt_f+engine.combustor.dmdt_f)*engine.nozzle_op.V^2/2-ambient.dmdt*ambient.V^2/2;
];
eta_th = DKE./[-engine.combustor.dmdt_f*1e3*engine.combustor.Q_R; -engine.combustor.dmdt_f*1e3*engine.combustor.Q_R-engine.afterburner_op.dmdt_f*engine.afterburner_op.Q_R*1e3];
eta_p = thrust.*ambient.V./DKE; % Propulsion efficiency, []
eta_0 = eta_th.*eta_p; % Overall efficiency, []
% eta_0 = thrust.*[engine.nozzle.V; engine.nozzle_op.V]./[-engine.combustor.dmdt_f*1e3*engine.combustor.Q_R; -engine.combustor.dmdt_f*1e3*engine.combustor.Q_R-engine.afterburner_op.dmdt_f*engine.afterburner_op.Q_R*1e3];

g=9.81;
LD = 7; % Estimate of cruise lift to drag ratio, []
s = -eta_0.*LD*engine.combustor.Q_R/g*log(1.66); % Aircraft range (assuming fuel is 40% of aircraft weight), km
% s = LD*ambient.V/g/TSFC*log(1.66)/1000;
table(thrust./1000, ST, TSFC, eta_th, eta_p, eta_0 , s,'VariableNames',{'Thrust [kN]','ST','TSFC','eta_th','eta_p','eta_0','s [km]'},'RowNames',{'No ABR','ABR'})

%}




R = 287;
ISA = load('ISA.mat');
ISA = ISA.ISA;

% M = 0.1:0.2:0.7;
% h = 0:5000:15000;
% ISADEV = -10:10:10;
% BPR = 0:0.5:1.5;
% HPCPR = 10:2:16;
% LPCPR = 1.2:0.5:2.2;
% TIT = 1750:250:2250;


M = 0.5;
h = 0;
ISADEV = 0;
BPR = 0.57;
HPCPR = 12.8;
LPCPR = 1.25;
TIT = 2000;

N = length(M)*length(h)*length(ISADEV)*length(BPR)*length(HPCPR)*length(LPCPR)*length(TIT);

for i = 1:length(M)
    for j = 1:length(h)
        for k = 1:length(ISADEV)
            for l = 1:length(BPR)
                for m = 1:length(HPCPR)
                    for n = 1:length(LPCPR)
                        for o = 1:length(TIT)
disp(['i=',num2str(i),' j=',num2str(j),' k=',num2str(k),' l=',num2str(l),' m=',num2str(m),' n=',num2str(n),' o=',num2str(o)])
% Ambient Conditions
ambient.T = interp1(ISA(:,1),ISA(:,2),h(j))+ISADEV(k);
ambient.p = interp1(ISA(:,1),ISA(:,3),h(j))*1e5;
ambient.rho = py.CoolProp.CoolProp.PropsSI('D','T',ambient.T,'P',ambient.p,'Air');
ambient.M = M(i);
cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',ambient.T,'P',ambient.p,'Air');
cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',ambient.T,'P',ambient.p,'Air');
ambient.h = py.CoolProp.CoolProp.PropsSI('H','T',ambient.T,'P',ambient.p,'Air');
ambient.gamma = cp/cv;
ambient.V = ambient.M*sqrt(ambient.gamma*R*ambient.T);
ambient.h0 = ambient.h + ambient.V^2/2;
ambient.p0 = ambient.p*(1+(ambient.gamma-1)/2*ambient.M^2)^(ambient.gamma/(ambient.gamma-1));
ambient.T0 = ambient.T*(1+(ambient.gamma-1)/2*ambient.M^2);

% Dry Engine Calculations
ambient.dmdt = 150;
dmdt_aH = ambient.dmdt/(BPR(l)+1);
dmdt_aC = BPR(l)*ambient.dmdt/(BPR(l)+1);
engine.dry.diffuser = diffuser(ambient, spr.INT, 0.95, R);
engine.dry.fan = comp(engine.dry.diffuser, FAN_pr, eta.FAN, 1.02, R);
dry.core = engine.dry.fan;
dry.core.dmdt = dmdt_aH;
dry.core = duct(dry.core, spr.LPC, R);
dry.bypass = engine.dry.fan;
dry.bypass.dmdt = dmdt_aC;
dry.bypass = duct(dry.bypass, spr.BPD, R);
engine.dry.lpc = comp(dry.core, LPCPR(n), eta.LPC, 3, R);
engine.dry.lpc_ducted = duct(engine.dry.lpc, spr.HPC, R);
engine.dry.hpc = comp(engine.dry.lpc_ducted, HPCPR(m), eta.HPC, 3, R); % Adjusted to get 150m/s in combustor in static conditions
engine.dry.combustor = combustor(engine.dry.hpc, spr.BRN, LHV, eta.BRN, R, TIT(o));
engine.dry.hpt = turb(engine.dry.combustor, engine.dry.hpc.w/eta.SFT, eta.HPT, 0.35, R);
engine.dry.lpt = turb(engine.dry.hpt, (engine.dry.lpc.w*dmdt_aH+engine.dry.fan.w*ambient.dmdt)/dmdt_aH/eta.SFT, eta.LPT, 0.7, R);
engine.dry.mixer = mix(engine.dry.lpt, dry.bypass, spr.MXR, R);
engine.dry.afterburner = afterburner_inop(engine.dry.mixer, spr.ABR, R);
engine.dry.nozzle = nozzle(engine.dry.afterburner, eta.NOZ, pi*(0.78/2)^2, R, ambient.p);

% Wet Engine Calculations
ambient.dmdt = 165;
dmdt_aH = ambient.dmdt/(BPR(l)+1);
dmdt_aC = BPR(l)*ambient.dmdt/(BPR(l)+1);
engine.wet.diffuser = diffuser(ambient, spr.INT, 0.95, R);
engine.wet.fan = comp(engine.wet.diffuser, FAN_pr, eta.FAN, 1.02, R);
wet.core = engine.wet.fan;
wet.core.dmdt = dmdt_aH;
wet.core = duct(wet.core, spr.LPC, R);
wet.bypass = engine.wet.fan;
wet.bypass.dmdt = dmdt_aC;
wet.bypass = duct(wet.bypass, spr.BPD, R);
engine.wet.lpc = comp(wet.core, LPCPR(n), eta.LPC, 3, R);
engine.wet.lpc_ducted = duct(engine.wet.lpc, spr.HPC, R);
engine.wet.hpc = comp(engine.wet.lpc_ducted, HPCPR(m), eta.HPC, 3, R); % Adjusted to get 150m/s in combustor in static conditions
engine.wet.combustor = combustor(engine.wet.hpc, spr.BRN, LHV, eta.BRN, R, TIT(o));
engine.wet.hpt = turb(engine.wet.combustor, engine.wet.hpc.w/eta.SFT, eta.HPT, 0.35, R);
engine.wet.lpt = turb(engine.wet.hpt, (engine.wet.lpc.w*dmdt_aH+engine.wet.fan.w*ambient.dmdt)/dmdt_aH/eta.SFT, eta.LPT, 0.7, R);
engine.wet.mixer = mix(engine.wet.lpt, wet.bypass, spr.MXR, R);
engine.wet.afterburner = afterburner_op(engine.wet.mixer, spr.ABRON, LHV);
% throat: 0.92m, exit: 1.15m
engine.wet.nozzle = nozzle(engine.wet.afterburner, eta.NOZWET, pi*(0.92/2)^2, R, ambient.p);

% Difference between thrust (dry and wet) and published values

% T-s diagram
tv = [ambient.T, engine.dry.diffuser.T, engine.dry.fan.T, engine.dry.lpc.T, engine.dry.hpc.T, engine.dry.combustor.T, engine.dry.hpt.T, engine.dry.lpt.T, engine.dry.mixer.T, engine.dry.afterburner.T, engine.dry.nozzle.T];
hv = [ambient.h, engine.dry.diffuser.h, engine.dry.fan.h, engine.dry.lpc.h, engine.dry.hpc.h, engine.dry.combustor.h, engine.dry.hpt.h, engine.dry.lpt.h, engine.dry.mixer.h, engine.dry.afterburner.h, engine.dry.nozzle.h];
pv = [ambient.p, engine.dry.diffuser.p, engine.dry.fan.p, engine.dry.lpc.p, engine.dry.hpc.p, engine.dry.combustor.p, engine.dry.hpt.p, engine.dry.lpt.p, engine.dry.mixer.p, engine.dry.afterburner.p, engine.dry.nozzle.p];
for i = 1:length(tv)
    sv(i) = py.CoolProp.CoolProp.PropsSI('S','T',tv(i),'P',pv(i),'Air');
end
figure
plot(sv,tv,'r*-')
hold on
plot([sv(1:3),py.CoolProp.CoolProp.PropsSI('S','T',dry.bypass.T,'P',dry.bypass.p,'Air'), sv(9)],[tv(1:3),dry.bypass.T,tv(9)],'b-')
plot([sv(end),sv(1)],[tv(end),tv(1)],'b--')
text(sv,tv+10, {'AMB','DIF','FAN','LPC','HPC','BNR','HPT','LPT','MIX','ABR','NOZ'})
text(py.CoolProp.CoolProp.PropsSI('S','T',dry.bypass.T,'P',dry.bypass.p,'Air'),dry.bypass.T+10,'BPD')
ylabel('T [K]')
xlabel('s, J/kg-K')
                        end
                    end
                end
            end
        end
    end
end
