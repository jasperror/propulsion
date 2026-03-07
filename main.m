% ENGR477 Propulsion
% Project 1: F135 Engine
% Jasper Palmer & Jackson Gilbert

clear; close all; clc;


% Meeting Todos
% Run a bunch of CEA simulations
% What to do about afterburning condition -> Another burner

%% Constants

% Design-point Operating Conditions
dmdt_a = 150; % Total air mass flow rate, kg/s
pr = 28; % Overall pressure ratio, []
T_5 = 2000; % TIT, K
BPR = 0.57; % Bypass ratio, []
FAN_pr = 1.75; % Fan pressure ratio, []
LPC_pr = 1.25; % LPC pressure ratio, []
HPC_pr = 12.8; % HPC pressure ratio, []
LHV = 43150; % Fuel heating value, kJ/kg
T_a = 288.15; % Design point ambient temperature, K
p_a = 101.325e3; % Ambient pressure, Pa

% Flow & Geometric Assumptions
M_1 = 0.5; % Fan inlet mach number, []
Vbar_45 = 150; % Combustor average axial velocity, m/s
d_9 = 0.78; % Nozzle throat diameter, m
d_10 = 0.78; % Nozzle exit diameter, m
M_7 = 0.5; % Turbine exit mach number, []

% Component Efficiencies
eta.FAN = 0.89; % Fan efficiency, []
eta.LPC = 0.88; % LPC efficiency, []
eta.HPC = 0.86; % HPC efficiency, []
eta.BRN = 0.99; % Burner efficiency, []
eta.HPT = 0.89; % HPT efficiency, []
eta.LPT = 0.91; % LPT efficiency, []
eta.NOZ = 0.98; % Nozzle efficiency, []
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
function val = cp_call(prop, n1, v1, n2, v2)
    val = py.CoolProp.CoolProp.PropsSI(prop, n1, v1, n2, v2, 'Air');
end
function pst(label, st, mdot)
    fprintf('%-22s  %8.1f   %8.2f   %9.3f\n', label, st.T0, st.p0/1e3, mdot);
end
CP = @(prop,name1,val1,name2,val2)py.CoolProp.CoolProp.PropsSI(prop,name1,val1,name2,val2,'Air');

function out = diffuser(in, spr_v, Ar, R)

    out.p0 = in.p0*spr_v;  out.T0 = in.T0;  out.dmdt = in.dmdt;
    T = in.T;  p = in.p;
    for k = 1:10000
        cp  = cp_call('CPMASS','T',T,'P',p);
        cv  = cp_call('CVMASS','T',T,'P',p);
        gam = cp/cv;
        rho = cp_call('D','T',T,'P',p);
        V   = in.V * in.rho / (rho * Ar);
        M   = V / sqrt(gam*R*T);
        pn  = out.p0*(1+(gam-1)/2*M^2)^(gam/(1-gam));
        Tn  = out.T0/(1+(gam-1)/2*M^2);
        err = sqrt(((Tn-T)/Tn)^2 + ((pn-p)/pn)^2);
        T = Tn;  p = pn;
        if err < 1e-7, break; end
    end
    out.T     = T;  out.p = p;
    out.rho   = cp_call('D','T',T,'P',p);
    cp = cp_call('CPMASS','T',T,'P',p);
    cv = cp_call('CVMASS','T',T,'P',p);
    out.gamma = cp/cv;
    out.V     = in.V * in.rho / (out.rho * Ar);
    out.M     = out.V / sqrt(out.gamma*R*T);
    out.h     = cp_call('H','T',T,'P',p);
    out.h0    = out.h + out.V^2/2;
end

function out = comp(in, pr, eta_c, Ar, R)
    out.p = in.p*pr;  out.dmdt = in.dmdt;
    h1    = cp_call('H','T',in.T,'P',in.p);
    h01   = h1 + in.V^2/2;
    s1    = cp_call('S','T',in.T,'P',in.p);
    T2    = in.T;
    for k = 1:10000
        rho  = cp_call('D','T',T2,'P',out.p);
        V2   = in.V * in.rho / (rho * Ar);
        h2s  = cp_call('H','P',out.p,'S',s1);
        h02s = h2s + V2^2/2;
        h02  = h01 + (h02s - h01)/eta_c;
        h2   = h02 - V2^2/2;
        T2n  = cp_call('T','H',h2,'P',out.p);
        err  = abs((T2n-T2)/T2n);
        T2   = T2n;
        if err < 1e-7, break; end
    end
    out.T  = T2;
    cp = cp_call('CPMASS','T',T2,'P',out.p);
    cv = cp_call('CVMASS','T',T2,'P',out.p);
    out.gamma = cp/cv;
    out.rho   = cp_call('D','T',T2,'P',out.p);
    out.V     = in.V * in.rho / (out.rho * Ar);
    out.h     = h2;
    out.h0    = h02;
    out.M     = out.V / sqrt(out.gamma*R*T2);
    out.p0    = out.p*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.T0    = T2*(1+(out.gamma-1)/2*out.M^2);
    out.w     = h02 - h01;  
end

function out = turb(in, w_req, eta_t, Ar, R)
    out.dmdt = in.dmdt;
    h1   = cp_call('H','T',in.T,'P',in.p);
    h01  = h1 + in.V^2/2;
    s1   = cp_call('S','T',in.T,'P',in.p);
    h02  = h01 - w_req;
    h02s = h01 - w_req/eta_t;

    T2 = in.T;  p2 = in.p;
    for k = 1:10000
        rho2  = cp_call('D','T',T2,'P',p2);
        V2    = in.V * in.rho / (rho2 * Ar);
        h2s   = h02s - V2^2/2;
        h2    = h02  - V2^2/2;
        p2n   = cp_call('P','H',h2s,'S',s1);
        T2n   = cp_call('T','H',h2,'P',p2n);
        err   = sqrt(((p2n-p2)/p2n)^2 + ((T2n-T2)/T2n)^2);
        p2 = p2n;  T2 = T2n;
        if err < 1e-7, break; end
    end
    out.T  = T2;  out.p = p2;
    cp = cp_call('CPMASS','T',T2,'P',p2);
    cv = cp_call('CVMASS','T',T2,'P',p2);
    out.gamma = cp/cv;
    out.rho   = cp_call('D','T',T2,'P',p2);
    out.V     = in.V * in.rho / (out.rho * Ar);
    out.M     = out.V / sqrt(out.gamma*R*T2);
    out.p0    = p2*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.T0    = T2*(1+(out.gamma-1)/2*out.M^2);
    out.h     = h2;
    out.h0    = h02;
end

function out = duct(in, spr_v, R)
    out.p0 = in.p0*spr_v;  out.T0 = in.T0;  out.dmdt = in.dmdt;
    T = in.T;  p = in.p;
    for k = 1:10000
        cp  = cp_call('CPMASS','T',T,'P',p);
        cv  = cp_call('CVMASS','T',T,'P',p);
        gam = cp/cv;
        rho = cp_call('D','T',T,'P',p);
        V   = in.V * in.rho / rho;          
        M   = V / sqrt(gam*R*T);
        pn  = out.p0*(1+(gam-1)/2*M^2)^(gam/(1-gam));
        Tn  = out.T0/(1+(gam-1)/2*M^2);
        err = sqrt(((Tn-T)/Tn)^2 + ((pn-p)/pn)^2);
        T = Tn;  p = pn;
        if err < 1e-7, break; end
    end
    out.T     = T;  out.p = p;
    out.rho   = cp_call('D','T',T,'P',p);
    cp = cp_call('CPMASS','T',T,'P',p);
    cv = cp_call('CVMASS','T',T,'P',p);
    out.gamma = cp/cv;
    out.V     = in.V * in.rho / out.rho;
    out.M     = out.V / sqrt(out.gamma*R*T);
    out.h     = cp_call('H','T',T,'P',p);
    out.h0    = out.h + out.V^2/2;
end

% %function out = CEARUN(p, T, CEA, Tv)
%     % Import CEA data
%     sz = size(CEA, [1, 3]);
%     Tad = reshape(CEA(:,4,:),sz);
%     c = reshape(CEA(:,7,:),sz);
%     h_BNR = reshape(CEA(:,5,:),sz);
%     gamma = reshape(CEA(:,6,:),sz);
%     CO = reshape(CEA(:,8,:),sz);
%     CO2 = reshape(CEA(:,9,:),sz);
%     NO = reshape(CEA(:,10,:),sz);
%     OF = reshape(CEA(1,2,:),[1, size(CEA,3)]);
%     phi = reshape(CEA(1,1,:),[1, size(CEA,3)]);
%     pv = CEA(:,3,1)';
%     h_air = py.CoolProp.CoolProp.PropsSI('H','T',T,'P',p,'Air')/1000-py.CoolProp.CoolProp.PropsSI('H','T',298.15,'P',101325,'Air')/1000; % Air absolute enthalpy [kJ/kg]
%     h_fuel = -284117/151.9; % Jet-A absolute enthalpy [kJ/kg]
% 
%     % Lookup Values
%     if T < 200 % If outside of bounds then use min/max. There is probably a better way to do this
%         out.OF = CEA(1,2,1);
%         out.phi = CEA(1,1,1);
%         h_reac = (h_air*OF + h_fuel)/(1+OF);
%         if p*1e-5 < 0.1
%             out.T = CEA(1,4,1);
%             out.c = CEA(1, 7, 1);
%             out.h = CEA(1, 5, 1);
%             out.Q_R = out.h-h_reac;
%             out.gamma = CEA(1,6,1);
%             out.CO2e = 10*CEA(1,8,1)+CEA(1,9,1);
%             out.NO = CEA(1,10,1);
%         elseif p*1e-5 > 50
%             out.T = CEA(end,4,1);
%             out.c = CEA(end, 7, 1);
%             out.h = CEA(end, 5, 1);
%             out.Q_R = out.h-h_reac;
%             out.gamma = CEA(end,6,1);
%             out.CO2e = 10*CEA(end,8,1)+CEA(end,9,1);
%             out.NO = CEA(end,10,1);
%         else
%             out.T = interp1(pv, CEA(:,4,1), p*1e-5);
%             out.c = interp1(pv, CEA(:,7,1), p*1e-5);
%             out.h = interp1(pv, CEA(:,5,1), p*1e-5);
%             out.Q_R = out.h-h_reac;
%             out.gamma = interp1(pv, CEA(:,6,1), p*1e-5);
%             out.CO2e = 10*interp1(pv, CEA(:,8,1), p*1e-5)+interp1(pv, CEA(:,9,1), p*1e-5);
%             out.NO = interp1(pv, CEA(:,10,1), p*1e-5);
%         end
%     elseif T > 2000
%         out.OF = CEA(1,2,end);
%         out.phi = CEA(1,1,end);
%         h_reac = (h_air*OF + h_fuel)/(1+OF);
%         if p*1e-5 < 0.1
%             out.T = CEA(1,4,end);
%             out.c = CEA(1, 7, end);
%             out.h = CEA(1, 5, end);
%             out.Q_R = out.h-h_reac;
%             out.gamma = CEA(1,6,end);
%             out.CO2e = 10*CEA(1,8,end)+CEA(1,9,end);
%             out.NO = CEA(1,10,end);
%         elseif p*1e-5 > 50
%             out.T = CEA(end,4,end);
%             out.c = CEA(end, 7, end);
%             out.h = CEA(end, 5, end);
%             out.Q_R = out.h-h_reac;
%             out.gamma = CEA(end,6,end);
%             out.CO2e = 10*CEA(end,8,end)+CEA(end,9,end);
%             out.NO = CEA(end,10,end);
%         else
%             out.T = interp1(pv, CEA(:, 4, end), p*1e-5);
%             out.c = interp1(pv, CEA(:, 7, end), p*1e-5);
%             out.h = interp1(pv, CEA(:, 5, end), p*1e-5);
%             out.Q_R = out.h-h_reac;
%             out.gamma = interp1(pv, CEA(:,6,end), p*1e-5);
%             out.CO2e = 10*interp1(pv, CEA(:,8,end), p*1e-5)+interp1(pv, CEA(:,9,end), p*1e-5);
%             out.NO = interp1(pv, CEA(:,10,end), p*1e-5);
%         end
%     else
%         out.OF = interp1(Tv, OF, T);
%         out.phi = interp1(Tv, phi, T);
%         h_reac = (h_air*OF + h_fuel)/(1+OF);
%         if p*1e-5 < 0.1
%             out.T = interp1(Tv, CEA(1,4,:), T);
%             out.c = interp1(Tv, CEA(1,7,:), T);
%             out.h = interp1(Tv, CEA(1,5,:), T);
%             out.Q_R = out.h-h_reac;
%             out.gamma = interp1(Tv, CEA(1,6,:), T);
%             out.CO2e = 10*interp1(Tv, CEA(1,8,:), T)+interp1(Tv, CEA(1,9,:), T);
%             out.NO = interp1(Tv, CEA(1,10,:), T);
%         elseif p*1e-5 > 50
%             out.T = interp1(Tv, CEA(end,4,:), T);
%             out.c = interp1(Tv, CEA(end,7,:), T);
%             out.h = interp1(Tv, CEA(end,5,:), T);
%             out.Q_R = out.h-h_reac;
%             out.gamma = interp1(Tv, CEA(end,6,:), T);
%             out.CO2e = 10*interp1(Tv, CEA(end,8,:), T)+interp1(Tv, CEA(end,9,:), T);
%             out.NO = interp1(Tv, CEA(end,10,:), T);
%         else % If parameters are OK then interpolate values from CEA
%             [X, Y] = meshgrid(Tv, pv);
%             out.T = interp2(X, Y, Tad, T, p*1e-5);
%             out.c = interp2(X, Y, c, T, p*1e-5);
%             out.h = interp2(X, Y, h_BNR, T, p*1e-5);
%             out.Q_R = out.h-h_reac;
%             out.gamma = interp2(X, Y, gamma, T, p*1e-5);
%             out.CO2e = 10*interp2(X, Y, CO, T, p*1e-5) + interp2(X, Y, CO2, T, p*1e-5);
%             out.NO = interp2(X, Y, NO, T, p*1e-5);
%         end
%     end
%     out.Q_R = -(py.CoolProp.CoolProp.PropsSI('H','T',out.T,'P',p,'Air')/1000-py.CoolProp.CoolProp.PropsSI('H','T',T,'P',p,'Air')/1000); % super scuffed not using the CEA (impove me!)
% end

function out = combustor(in, spr_v, LHV, eta_b, R)

    TIT    = 2000;   % K
    out    = in;
    out.p0 = spr_v * in.p0;
    out.T0 = TIT;

    h04 = in.h0;
    h05 = cp_call('H','T',TIT,'P',out.p0) + in.V^2/2;

    out.dmdt_f = in.dmdt * (h05 - h04) / (eta_b * LHV * 1000);
    out.dmdt   = in.dmdt + out.dmdt_f;
    out.V      = in.V;
    out.h0     = h05;
    out.Q_R    = LHV;  % kJ/kg_fuel, stored for reference

    cp = cp_call('CPMASS','T',TIT,'P',out.p0);
    cv = cp_call('CVMASS','T',TIT,'P',out.p0);
    out.gamma = cp/cv;
    out.M     = out.V / sqrt(out.gamma*R*TIT);
    out.T     = TIT / (1+(out.gamma-1)/2*out.M^2);
    out.p     = out.p0*(1+(out.gamma-1)/2*out.M^2)^(-out.gamma/(out.gamma-1));
    out.rho   = cp_call('D','T',out.T,'P',out.p);
    out.h     = cp_call('H','T',out.T,'P',out.p);
end


function out = mix(core, byp, spr_v, R)
    out.dmdt = core.dmdt + byp.dmdt;
    out.h0   = (core.dmdt*core.h0 + byp.dmdt*byp.h0) / out.dmdt;
    out.p    = (core.p + byp.p) / 2;
    out.V    = (core.dmdt*core.V + byp.dmdt*byp.V) / out.dmdt;

    h = out.h0 - out.V^2/2;
    T = (core.T + byp.T) / 2;
    for k = 1:100
        h_k  = cp_call('H','T',T,'P',out.p);
        cp_k = cp_call('CPMASS','T',T,'P',out.p);
        dT   = (h - h_k) / cp_k;
        T    = T + dT;
        if abs(dT) < 0.005, break; end
    end
    out.T   = T;
    out.rho = cp_call('D','T',T,'P',out.p);
    cp = cp_call('CPMASS','T',T,'P',out.p);
    cv = cp_call('CVMASS','T',T,'P',out.p);
    out.gamma = cp/cv;
    out.M     = out.V / sqrt(out.gamma*R*T);
    out.T0    = T*(1+(out.gamma-1)/2*out.M^2);
    out.p0    = spr_v*out.p*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.h     = h;
end

function out = afterburner_inop(in, spr_v, R)
    out        = duct(in, spr_v, R);
    out.dmdt_f = 0;
end

%% ---- AFTERBURNER (on) ---------------------------------------------------
function out = afterburner_op(in, spr_v, LHV, eta_b, R)
    T09    = 2450;   % K
    out    = in;
    out.p0 = spr_v * in.p0;

    h08 = in.h0;
    h09 = cp_call('H','T',T09,'P',out.p0) + in.V^2/2;

    out.dmdt_f = in.dmdt * (h09 - h08) / (eta_b * LHV * 1000);
    out.dmdt   = in.dmdt + out.dmdt_f;
    out.T0     = T09;
    out.h0     = h09;
    out.V      = in.V;

    cp = cp_call('CPMASS','T',T09,'P',out.p0);
    cv = cp_call('CVMASS','T',T09,'P',out.p0);
    out.gamma = cp/cv;
    out.M     = out.V / sqrt(out.gamma*R*T09);
    out.T     = T09 / (1+(out.gamma-1)/2*out.M^2);
    out.p     = out.p0*(1+(out.gamma-1)/2*out.M^2)^(-out.gamma/(out.gamma-1));
    out.rho   = cp_call('D','T',out.T,'P',out.p);
    out.h     = cp_call('H','T',out.T,'P',out.p);
end

%% ---- NOZZLE -------------------------------------------------------------
function out = nozzle(in, eta_n, A, R, p_amb)

    out   = in;
    mdot  = in.dmdt;
    h1    = cp_call('H','T',in.T,'P',in.p);
    h01   = h1 + in.V^2/2;
    gam   = in.gamma;
    choke = ((gam+1)/2)^(gam/(gam-1));

    if in.p0/p_amb >= choke
        %--- Choked ---
        out.M = 1;
        out.T = in.T0 / (1+(gam-1)/2);
        out.p = in.p0 * (2/(gam+1))^(gam/(gam-1));
        cp = cp_call('CPMASS','T',out.T,'P',out.p);
        cv = cp_call('CVMASS','T',out.T,'P',out.p);
        out.gamma = cp/cv;
        a     = sqrt(out.gamma*R*out.T);
        out.V = a;
        out.rho = cp_call('D','T',out.T,'P',out.p);
       
        out.h  = cp_call('H','T',out.T,'P',out.p);
        out.h0 = out.h + a^2/2;
    else
        %--- Unchoked ---
        out.p = p_amb;
        Tg    = in.T * 0.9;
        for k = 1:1000
            rho  = cp_call('D','T',Tg,'P',out.p);
            V    = mdot / (rho*A);
            h2   = h01 - eta_n*V^2/2;
            Tnew = cp_call('T','H',h2,'P',out.p);
            err  = abs((Tnew-Tg)/Tnew);
            Tg   = Tnew;
            if err < 1e-8, break; end
        end
        out.T   = Tnew;  out.V = V;
        out.rho = cp_call('D','T',Tnew,'P',out.p);
        cp = cp_call('CPMASS','T',Tnew,'P',out.p);
        cv = cp_call('CVMASS','T',Tnew,'P',out.p);
        out.gamma = cp/cv;
        out.M     = V / sqrt(out.gamma*R*Tnew);
        out.h     = h2;
    end
    out.p0 = out.p*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.T0 = out.T*(1+(out.gamma-1)/2*out.M^2);
end

    % i = 0; err = 1; tmp.T = in.T0; tmp.M = 0.5;
    % disp(['p0=',num2str(out.p0*1e-5),', T0=',num2str(out.T0)])
    % while err > 0.001 && i < 1e4
    %     i = i+1;
    %     disp(['i=',num2str(i),', T=',num2str(tmp),', p=',num2str(p)])
    %     [~, ~, ~, gamma] = air_props(tmp.T);
    %     out.p = out.p0*(1+(gamma-1)/2*tmp.M^2)^(-gamma/(gamma-1));
    %     out.u = sqrt(2*eta*gamma/(gamma-1)*R*in.T0*(1-(p/in.p0)^((gamma-1)/gamma)));
    %     M = out.u/sqrt(gamma*R*tmp.T);
    %     T = in.T0/(1+(gamma-1)/2*M^2);
    %     err = sqrt(mean( ((T-tmp.T)/T)^2+((M-tmp.M)/M)^2 ));
    %     tmp.T = T;
    %     tmp.M = M;
    % end



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

% Ambient Conditions
ambient.T = 15 + 273.15;   % K
ambient.p = 101.325e3;     % Pa
ambient.dmdt = 150;           % kg/s total air
ambient.M  = 0.5;
R = 287;

ambient.rho = CP('D','T',ambient.T,'P',ambient.p);
cp_a = CP('CPMASS','T',ambient.T,'P',ambient.p);
cv_a = CP('CVMASS','T',ambient.T,'P',ambient.p);
ambient.gamma = cp_a / cv_a;
ambient.V = ambient.M * sqrt(ambient.gamma * R * ambient.T);
ambient.h = CP('H','T',ambient.T,'P',ambient.p);
ambient.h0 = ambient.h + ambient.V^2/2;
ambient.p0 = ambient.p * (1+(ambient.gamma-1)/2*ambient.M^2)^(ambient.gamma/(ambient.gamma-1));
ambient.T0 = ambient.T * (1+(ambient.gamma-1)/2*ambient.M^2);

dmdt_aH = ambient.dmdt / (BPR+1);        % core mass flow, kg/s
dmdt_aC = BPR * ambient.dmdt / (BPR+1);  % bypass mass flow, kg/s

% Engine Calculations
engine.diffuser = diffuser(ambient, spr.INT, 0.95, R);
engine.fan      = comp(engine.diffuser, FAN_pr, eta.FAN, 1.02, R);

% Fan exit split
core        = engine.fan;  core.dmdt   = dmdt_aH;
bypass      = engine.fan;  bypass.dmdt = dmdt_aC;
core        = duct(core,   spr.LPC, R);
bypass      = duct(bypass, spr.BPD, R);

% Compressors
engine.lpc = comp(core, LPC_pr, eta.LPC, 3, R);
engine.lpc_ducted = duct(engine.lpc, spr.HPC, R);
engine.hpc = comp(engine.lpc_ducted, HPC_pr, eta.HPC, 3, R);

% Combustor (targets TIT = 2000 K)
engine.combustor = combustor(engine.hpc, spr.BRN, LHV, eta.BRN, R);

% HP shaft power balance
w_HPT = engine.hpc.w * engine.lpc_ducted.dmdt / (eta.SFT * engine.combustor.dmdt);
engine.hpt = turb(engine.combustor, w_HPT, eta.HPT, 0.35, R);

% LP shaft power balance

P_LP = engine.fan.w * ambient.dmdt + engine.lpc.w * dmdt_aH;
w_LPT = P_LP / (eta.SFT * engine.hpt.dmdt);
engine.lpt = turb(engine.hpt, w_LPT, eta.LPT, 0.7, R);

% Mixer
engine.mixer = mix(engine.lpt, bypass, spr.MXR, R);

% Afterburner
engine.afterburner = afterburner_inop(engine.mixer, spr.ABR, R);
engine.afterburner.Q_R = 0;
engine.afterburner_op = afterburner_op(engine.mixer, spr.ABRON, LHV, eta.BRN, R);

% Nozzles
A_dry = pi * (d_9/2)^2;
A_wet = pi * (0.92/2)^2;

engine.nozzle = nozzle(engine.afterburner, eta.NOZ, A_dry, R, ambient.p);
engine.nozzle_op = nozzle(engine.afterburner_op, 0.97, A_wet, R, ambient.p);

%output paramegters
g = 9.81;

% Mass flows at nozzle exit
mdot_e_dry = engine.afterburner.dmdt;
mdot_e_wet = engine.afterburner_op.dmdt;
mdot_f_dry = engine.combustor.dmdt_f;
mdot_f_wet = engine.combustor.dmdt_f + engine.afterburner_op.dmdt_f;

% Thrust [N]
thrust = [
    mdot_e_dry * engine.nozzle.V    - ambient.dmdt*ambient.V + A_dry*(engine.nozzle.p    - ambient.p);
    mdot_e_wet * engine.nozzle_op.V - ambient.dmdt*ambient.V + A_wet*(engine.nozzle_op.p - ambient.p)
];

% Specific thrust [N/(kg/s)]
ST = thrust ./ ambient.dmdt;

% TSFC [kg/(N.s)]
TSFC = [mdot_f_dry; mdot_f_wet] ./ thrust;

% Change in kinetic energy rate [W]
DKE = [
    0.5*mdot_e_dry*engine.nozzle.V^2 - 0.5*ambient.dmdt*ambient.V^2;
    0.5*mdot_e_wet*engine.nozzle_op.V^2 - 0.5*ambient.dmdt*ambient.V^2
];

% Thermal efficiency
eta_th = DKE ./ ([mdot_f_dry; mdot_f_wet] * LHV * 1000);

% Propulsive efficiency
eta_p = thrust .* ambient.V ./ DKE;

% Overall efficiency
eta_0 = eta_th .* eta_p;

% Breguet range [km]
LD = 7;
s  = ambient.V/g ./ TSFC * LD * log(1.66) / 1000;

% Station T0 / p0 summary
fprintf('\n===== STATION SUMMARY ============================\n')
fprintf('%-22s  %8s   %8s   %9s\n','Station','T0 [K]','p0 [kPa]','mdot [kg/s]')
pst('Ambient',        ambient,              ambient.dmdt)
pst('Fan exit',       engine.fan,           ambient.dmdt)
pst('LPC exit',       engine.lpc,           dmdt_aH)
pst('HPC exit',       engine.hpc,           dmdt_aH)
pst('Combustor',      engine.combustor,     engine.combustor.dmdt)
pst('HPT exit',       engine.hpt,           engine.hpt.dmdt)
pst('LPT exit',       engine.lpt,           engine.lpt.dmdt)
pst('Bypass duct',    bypass,               dmdt_aC)
pst('Mixer exit',     engine.mixer,         engine.mixer.dmdt)
pst('AB exit (dry)',  engine.afterburner,   engine.afterburner.dmdt)
pst('AB exit (wet)',  engine.afterburner_op,engine.afterburner_op.dmdt)
pst('Nozzle (dry)',   engine.nozzle,        mdot_e_dry)
pst('Nozzle (wet)',   engine.nozzle_op,     mdot_e_wet)

% Performance table
fprintf('\n===== PERFORMANCE SUMMARY =======================\n')
disp(table(thrust./1000, ST, TSFC, eta_th, eta_p, eta_0, s, ...
    'VariableNames',{'Thrust_kN','ST','TSFC','eta_th','eta_p','eta_0','Range_km'}, ...
    'RowNames',{'Dry','Wet_AB'}))

fprintf('Fuel flow dry  : %.4f kg/s\n', mdot_f_dry)
fprintf('Fuel flow wet  : %.4f kg/s\n', mdot_f_wet)
fprintf('Nozzle V dry   : %.1f m/s   p = %.2f kPa\n', engine.nozzle.V,    engine.nozzle.p/1e3)
fprintf('Nozzle V wet   : %.1f m/s   p = %.2f kPa\n', engine.nozzle_op.V, engine.nozzle_op.p/1e3)
