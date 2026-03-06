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
LPC_pr = 1.25; % LPC pressure ratio, []
T_5 = 2000; % TIT, K
BPR = 0.57; % Bypass ratio, []
FAN_pr = 1.75; % Fan pressure ratio, []
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
spr.NOZ = 0.98; % Nozzle SPR, []

%% Inlet ambient conditions

% gamma_guess = 1.4;
% flow.p0 = p_a * (1 + ((gamma_guess - 1)/2) * M_1^2)^(gamma_guess/(gamma_guess - 1)); % Ambiant static pressure, Pa
% flow.T0 = T_a * (1 + ((gamma_guess - 1)/2) * M_1^2); % Ambiant static temperature, K
% flow.mdot = dmdt_a; % Mass flow rate, kg/s
% flow.M = M_1;
% flow.T = T_a;
% flow.p = p_a;
% flow.V = M_1*sqrt(gamma_guess*R*T_a);
% [flow.cp, flow.h, flow.s] = air_props(flow.T0);

%% Engine Components

function [cp,h,s,gamma] = air_props(T)
    R = 0.287;  % KJ/kg-K
    
    if T <= 1000
        a1 = 3.78245636;
        a2 = -2.99673416e-3;
        a3 = 9.84730201e-6;
        a4 = -9.68129509e-9;
        a5 = 3.24372837e-12;
        a6 = -1063.94356;
        a7 = 3.65767573;
    else
        a1 = 3.28253784;
        a2 = 1.48308754e-3;
        a3 = -7.57966669e-7;
        a4 = 2.09470555e-10;
        a5 = -2.16717794e-14;
        a6 = -1088.45772;
        a7 = 5.45323129;
    end
    
    cp = R*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);
    h  = R*T*(a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T);
    s  = R*(a1*log(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7);
    cv = cp - R;
    gamma = cp/cv;
end

function s = entropy_air(T)
    [~,~,s,~] = air_props(T);
end

function h = enthalpy_air(T)
    [~,h,~,~] = air_props(T);
end

function out = diffuser(in, spr, Ar) % Assuming adiabatic
    out = in;
    out.p0 = in.p0 * spr;
    [~, ~, ~, gamma1] = air_props(in.T);

    i = 0; err = 1; tmp = in.T0;
    while err > 0.001 && i < 1e4
        % disp(['DIFFUSER i=',num2str(i),', T=',num2str(tmp)])
        i = i+1;
        [~, ~, ~, gamma2] = air_props(tmp);
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',tmp,'P',out.p,'Air');
        out.M = Ar*in.rho/out.rho*in.M*sqrt(gamma1*in.T)/sqrt(gamma2*tmp); % From mass conservation
        out.p = out.p0*(1+(gamma2-1)/2*out.M^2)^(gamma2/(1-gamma2));
        out.T = in.T0/(1+(gamma2-1)/2*out.M^2);
        err = abs((out.T-tmp)/out.T);
        tmp = out.T;
    end
    out.V = in.rho/out.rho*Ar*in.V;
    [~, ~, ~, out.gamma] = air_props(out.T);
end

function out = comp(in, pr, eta, Ar, R)
    out.p = in.p*pr;
    out.dmdt = in.dmdt;
    
    tmp = in.T; i = 0; err = 1;
    while err > 0.001 && i < 1e4
        i = i + 1;
        
        h1 = py.CoolProp.CoolProp.PropsSI('H','T',in.T,'P',in.p,'Air');
        h01 = h1 + in.V^2/2;
        s1 = py.CoolProp.CoolProp.PropsSI('S','T',in.T,'P',in.p,'Air');
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
    out.M = out.V/sqrt(out.gamma*R*out.T);
    out.p0 = out.p*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.T0 = out.T*(1+(out.gamma-1)/2*out.M^2);
    out.w = h2-h1+0.5*(out.V^2-in.V^2);
end

function out = turb(in, pr, eta, Ar, R)
    out.p = in.p*pr;
    out.dmdt = in.dmdt;
    
    tmp = in.T; i = 0; err = 1;
    while err > 0.001 && i < 1e4

        i = i + 1;
        
        h1 = py.CoolProp.CoolProp.PropsSI('H','T',in.T,'P',in.p,'Air');
        h01 = h1 + in.V^2/2;
        s1 = py.CoolProp.CoolProp.PropsSI('S','T',in.T,'P',in.p,'Air');
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',tmp,'P',out.p,'Air');
        out.V = in.rho/out.rho*Ar*in.V;
        h2s = py.CoolProp.CoolProp.PropsSI('H','P',out.p,'S',s1,'Air');
        h02s = h2s + out.V^2/2;
        h02 = h01-(h01-h02s)*eta;
        h2 = h02 - out.V^2/2;
        out.T = py.CoolProp.CoolProp.PropsSI('T','H',h2,'P',out.p,'Air');

        err = abs((out.T-tmp)/out.T);
        tmp = out.T;
    end

    
    out.M = out.V/sqrt(out.gamma*R*out.T);
    out.p0 = out.p*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(out.gamma-1));
    out.T0 = out.T*(1+(out.gamma-1)/2*out.M^2);
    out.w = h2-h1+0.5*(out.V^2-in.V^2);
end

function out = duct(in, spr)
    out.dmdt = in.dmdt;
    out.p0 = in.p0 * spr;
    out.T0 = in.T0;
    tmp = in.T; i = 0; err = 1;
    while err > 0.001 && i < 1e4
        i = i + 1;
        cp = py.CoolProp.CoolProp.PropsSI('CPMASS','T',tmp,'P',out.p,'Air');
        cv = py.CoolProp.CoolProp.PropsSI('CVMASS','T',tmp,'P',out.p,'Air');
        out.gamma = cp/cv;
        out.rho = py.CoolProp.CoolProp.PropsSI('D','T',tmp,'P',out.p,'Air');
        out.V = in.rho/out.rho*in.V; % Assuming no change in cross-section area
        out.M = in.V/sqrt(out.gamma*R*tmp);
        out.p = out.p0*(1+(out.gamma-1)/2*out.M^2)^(out.gamma/(1-out.gamma));
        out.T = out.T0/(1+(gamma-1)/2*out.M^2);
        err = abs((out.T-tmp)/out.T);
        tmp = out.T;
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
    pv = CEA(:,3,1)';
    [~, h_air] = air_props(T);
    h_fuel = -284117/151.9; % Jet-A enthalpy [kJ/kg]

    % Lookup Values
    if T < 200 % If outside of bounds then use min/max. There is probably a better way to do this
        out.OF = CEA(1,2,1);
        h_reac = (h_air*OF + h_fuel)/(1+OF);
        if p*1e-5 < 0.1
            out.T = CEA(1,4,1);
            out.c = CEA(1, 7, 1);
            out.Q_R = CEA(1, 5, 1)-h_reac;
            out.gamma = CEA(1,6,1);
            out.CO2e = 10*CEA(1,8,1)+CEA(1,9,1);
            out.NO = CEA(1,10,1);
        elseif p*1e-5 > 50
            out.T = CEA(end,4,1);
            out.c = CEA(end, 7, 1);
            out.Q_R = CEA(end, 5, 1)-h_reac;
            out.gamma = CEA(end,6,1);
            out.CO2e = 10*CEA(end,8,1)+CEA(end,9,1);
            out.NO = CEA(end,10,1);
        else
            out.T = interp1(pv, CEA(:,4,1), p*1e-5);
            out.c = interp1(pv, CEA(:,7,1), p*1e-5);
            out.Q_R = interp1(pv, CEA(:,5,1), p*1e-5)-h_reac;
            out.gamma = interp1(pv, CEA(:,6,1), p*1e-5);
            out.CO2e = 10*interp1(pv, CEA(:,8,1), p*1e-5)+interp1(pv, CEA(:,9,1), p*1e-5);
            out.NO = interp1(pv, CEA(:,10,1), p*1e-5);
        end
    elseif T > 2000
        out.OF = CEA(1,2,end);
        h_reac = (h_air*OF + h_fuel)/(1+OF);
        if p*1e-5 < 0.1
            out.T = CEA(1,4,end);
            out.c = CEA(1, 7, end);
            out.Q_R = CEA(1, 5, end)-h_reac;
            out.gamma = CEA(1,6,end);
            out.CO2e = 10*CEA(1,8,end)+CEA(1,9,end);
            out.NO = CEA(1,10,end);
        elseif p*1e-5 > 50
            out.T = CEA(end,4,end);
            out.c = CEA(end, 7, end);
            out.Q_R = CEA(end, 5, end)-h_reac;
            out.gamma = CEA(end,6,end);
            out.CO2e = 10*CEA(end,8,end)+CEA(end,9,end);
            out.NO = CEA(end,10,end);
        else
            out.T = interp1(pv, CEA(:, 4, end), p*1e-5);
            out.c = interp1(pv, CEA(:, 7, end), p*1e-5);
            out.Q_R = interp1(pv, CEA(:, 5, end), p*1e-5)-h_reac;
            out.gamma = interp1(pv, CEA(:,6,end), p*1e-5);
            out.CO2e = 10*interp1(pv, CEA(:,8,end), p*1e-5)+interp1(pv, CEA(:,9,end), p*1e-5);
            out.NO = interp1(pv, CEA(:,10,end), p*1e-5);
        end
    else
        out.OF = interp1(Tv, OF, T);
        h_reac = (h_air*OF + h_fuel)/(1+OF);
        if p*1e-5 < 0.1
            out.T = interp1(Tv, CEA(1,4,:), T);
            out.c = interp1(Tv, CEA(1,7,:), T);
            out.Q_R = interp1(Tv, CEA(1,5,:), T)-h_reac;
            out.gamma = interp1(Tv, CEA(1,6,:), T);
            out.CO2e = 10*interp1(Tv, CEA(1,8,:), T)+interp1(Tv, CEA(1,9,:), T);
            out.NO = interp1(Tv, CEA(1,10,:), T);
        elseif p*1e-5 > 50
            out.T = interp1(Tv, CEA(end,4,:), T);
            out.c = interp1(Tv, CEA(end,7,:), T);
            out.Q_R = interp1(Tv, CEA(end,5,:), T)-h_reac;
            out.gamma = interp1(Tv, CEA(end,6,:), T);
            out.CO2e = 10*interp1(Tv, CEA(end,8,:), T)+interp1(Tv, CEA(end,9,:), T);
            out.NO = interp1(Tv, CEA(end,10,:), T);
        else % If parameters are OK then interpolate values from CEA
            [X, Y] = meshgrid(Tv, pv);
            size(pv)
            size(Tv)
            size(X)
            out.T = interp2(X, Y, Tad, p*1e-5, T);
            out.c = interp2(X, Y, c, p*1e-5, T);
            out.Q_R = interp2(X, Y, h_BNR, p*1e-5, T)-h_reac; % There shouldn't really need to be an absolute value here but it is hard to get enthalpy for fuel
            out.gamma = interp2(X, Y, gamma, p*1e-5, T);
            out.CO2e = 10*interp2(X, Y, CO, p*1e-5, T) + interp2(X, Y, CO2, p*1e-5, T);
            out.NO = interp2(X, Y, NO, p*1e-5, T);
        end
    end
    out.Q_R = 1e3*out.Q_R;
end

function out = combustor(in, spr, LHV, eta, R)
    i = 0; err = 1; tmp = in.T0;
    while err > 0.001 && i < 1e4
        disp(['i=',num2str(i),', T=',num2str(tmp)])
        i = i+1;
        [~, ~, ~, gamma] = air_props(tmp);
        M = 250/sqrt(gamma*R*tmp);
        T = in.T0/(1+(gamma-1)/2*M^2);
        err = abs((T-tmp)/T);
        tmp = T;
    end
    p = in.p0*(1+(gamma-1)/2*M^2)^(-gamma/(gamma-1));

    CEA = load('CEA.mat');
    out = CEARUN(p, T, CEA.COMB, [200, 500,
        700, 800, 900, 1000, 2000]);
    out.p = p;
    out.p0 = spr*in.p0;

    % Calculate resulting parameters
    out.V = 250;
    out.M = out.V/out.c;
    t0_star = in.T0*(1+out.gamma*out.M^2)^2/(2*(out.gamma+1)*out.M^2*(1+(out.gamma-1)/2*out.M^2));
    out.T0 = t0_star*(in.T0/t0_star + out.Q_R/(out.gamma*R/(out.gamma-1)*t0_star));
    out.dmdt_f = out.Q_R/LHV/eta;
    out.dmdt = out.dmdt_f+in.dmdt;
end

function out = mixer(core, bypass, spr)
% Assuming 100% mixing takes place
  out = core;
  mdot_core = core.dmdt;
  mdot_bypass = bypass.dmdt;
  mdot_a = mdot_bypass + mdot_core;

  [~,hc,~,~] = air_props(core.T);
  [~,hb,~,~] = air_props(bypass.T);

  hm = (mdot_core * hc + mdot_bypass * hb) / mdot_a;

  fun = @(T) enthalpy_air(T) - hm;
  T_mix = fzero(fun, (core.T + bypass.T)/2);

  p_mix = min(core.p, bypass.p);
  out.p = p_mix * spr;
  out.T = T_mix;
  out.dmdt = mdot_a;
end

function out = afterburner_inop(in, spr)
    out = duct(in, spr);
    out.dmdt_f = 0;
end

function out = afterburner_op(in, spr, LHV, eta)
    CEA = load('CEA.mat');
    out = CEARUN(in.p, in.T, CEA.ABR, [200, 500, 700, 900, 2000]);
    out.p0 = in.p0*spr;
    out.dmdt_f = out.Q_R/LHV/eta;
end

function out = nozzle(in, spr, eta, R)
% need to check for choked flow
    out.h0 = in.h+in.V^2/2;
    out.h = in.h;
    out.p0 = in.p0*spr;
    out.T0 = in.T0;
    M = in.M;

    [~,~,~,gamma] = air_props(in.T);
    CPR = (2/(gamma+1))^(gamma/(gamma-1));
    actual_PR = (2/(gamma+1))^(gamma/(gamma-1)) * ((1+gamma)/(1+gamma*in.M^2)) * (1+(gamma-1)/2*in.M^2)^(gamma/(gamma-1));
    if CPR >= actual_PR
        warning('Flow is chocked')
        T0_star = T0 / ((2*(gamma+1)*M^2 * (1 + ((gamma-1)/2)*M^2)) / (1 + gamma*M^2)^2);
        p0_star = in.p0*CPR;
        T_star = in.T / (((1 + gamma) / (1 + gamma*M^2))^2 * M^2);
        p_star = in.p / ((1 + gamma) / (1 + gamma*M^2));
        out.T = T_star;
        out.T0 = T0_star;
        out.p = p_star;
        out.p0 = p0_star;
        out.M = 1;
        out.u = sqrt(gamma * R * T_star);
        return
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
Jackson's organization (for reference for now)
a = flow;

st1 = diffuser(a, spr.int);
st2 = fan(st1, FAN_pr, eta.fan);
st2 = duct(st2,spr.LPC);
st3 = compressor(st2, LPC_pr, eta.Lpc);
st3 = duct(st3, spr.LPC);
st4 = compressor(st3, HPC_pr, eta.HPC);
st4 = duct(st4, spr.HPC);
st5 = combuster(st4, spr.BRN, Vbar_45);
st5 = duct(st5,spr.BRN);

W_HPC = st4.h - st3.h;
W_LPC = st3.h - st2.h;
W_FAN = st2.h - st1.h;
W_shaft1 = (W_LPC + W_FAN) / eta.SFT;

st6 = turbine(st5, W_HPC, eta.HPT);
%}

ambient.T = 15+273.15;
ambient.p = 101.325e3;
ambient.dmdt = 150;
ambient.rho = py.CoolProp.CoolProp.PropsSI('D','T',ambient.T,'P',ambient.p,'Air');
R = 287;
[~, ~, ~, ambient.gamma] = air_props(ambient.T);
ambient.M = 0.5;
ambient.V = ambient.M*sqrt(ambient.gamma*R*ambient.T);
dmdt_aH = ambient.dmdt/(BPR+1);
dmdt_aC = BPR*ambient.dmdt/(BPR+1);
u = ambient.M*sqrt(ambient.gamma*R*ambient.T);
ambient.p0 = ambient.p*(1+(ambient.gamma-1)/2*ambient.M^2)^(ambient.gamma/(ambient.gamma-1));
ambient.T0 = ambient.T*(1+(ambient.gamma-1)/2*ambient.M^2);



engine.diffuser = diffuser(ambient, spr.INT, 0.95);
engine.fan = comp(engine.diffuser, FAN_pr, eta.FAN, 1.02, R)
% Should account for duct p0 losses here
engine.lpc = comp(engine.fan, LPC_pr, eta.LPC, 1.5, R)
% and here
engine.hpc = comp(engine.lpc, HPC_pr, eta.HPC, 2.3, R)
% and here
engine.combustor = combustor(engine.hpc, spr.BRN, LHV, eta.BRN, R)
% and here
engine.hpt = turb(engine.combustor, HPT_pr, eta.HPT, 0.35, R)
engine.lpt = turb(engine.hpt, LPT_pr, eta.LPT, 0.7, R)
% % and here
% engine.mixer = mixer()
% % and here
% engine.afterburner = afterburner_inop()
% % and here
% engine.nozzle = nozzle()

%{
engine.diffuser = diffuser(ambient, spr.INT, Ar);
engine.fan = fan(engine.diffuser, FAN_pr, eta.FAN);

core = engine.fan;
bypass = engine.fan;
core.dmdt = ambient.dmdt / (1 + BPR);
bypass.dmdt = ambient.dmdt - core.dmdt;
bypass = duct(bypass, spr.BPD); 
core = duct(core, spr.LPC);

engine.lpc = compressor(core, spr.LPC, eta.LPC);
engine.hpc = compressor(engine.lpc, spr.HPC ,eta.HPC);

engine.hpc.T0 = 800;

engine.combustor = combustor(engine.hpc, spr.BRN, LHV, eta.BRN, R);
engine.hpt = turbine(engine.combustor, engine.hpc.w/eta.SFT, eta.HPT);
engine.lpt = turbine(engine.hpt, (core.dmdt*engine.lpc.w+ambient.dmdt*engine.fan.w)/core.dmdt/eta.SFT, eta.LPT);
engine.mixer = mixer(engine.lpt, bypass, spr.MXR);
engine.afterburner = afterburner_inop(engine.mixer, spr.ABR);
engine.nozzle = nozzle(engine.afterburner, spr.NOZ, eta.NOZ, R);

% parameters needed
p_4 = engine.hpc.p;
T_4 = engine.hpc.T;
dmdt_f = engine.combustor.dmdt_f+engine.afterburner.dmdt_f;
Q_R = engine.combustor.Q_R;
u_e = engine.nozzle.u;
p_e = engine.nozzle.p;
% derived parameters
f = dmdt_f/ambient.dmdt;

% Output Parameters
thrust = (ambient.dmdt+dmdt_f)*u_e-ambient.dmdt*u+pi*d_10^2/4*(p_e-p_a); % Thrust, N
ST = thrust/ambient.dmdt; % Specific thrust, Ns/kg
TSFC = dmdt_f/thrust; % Thrust-specific fuel consumption, kg/Ns

eta_th = (1+f)*(u_e^2-u^2)/(2*f*Q_R); % Thermal efficiency, []
eta_p = thrust*u/ambient.dmdt/((1+f)*(u_e^2/2)-u^2/2); % Propulsion efficiency, []
eta_0 = eta_th*eta_p; % Overall efficiency, []

LD = 10; % Estimate of cruise lift to drag ratio, []
s = eta_0*LD*Q_R/g*log(1.66); % Aircraft range (assuming fuel is 40% of aircraft weight), m
%}