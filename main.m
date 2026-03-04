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

gamma_guess = 1.4;
flow.p0 = p_a * (1 + ((gamma_guess - 1)/2) * M_1^2)^(gamma_guess/(gamma_guess - 1)); % Ambiant static pressure, Pa
flow.T0 = T_a * (1 + ((gamma_guess - 1)/2) * M_1^2); % Ambiant static temperature, K
flow.mdot = dmdt_a; % Mass flow rate, kg/s
flow.M = M_1;
flow.T = T_a;
flow.p = p_a;
flow.V = M_1*sqrt(gamma_guess*R*T_a);
[flow.cp, flow.h, flow.s] = air_props(flow.T0);

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

function out = diffuser(in,pr) %Isentropic and adiabatic
    out = in;
    out.p0 = in.p0 * pr; % does this need to use spr or is pr fine
    out.T = in.T;
    out.p = in.p * pr;
     
end

function out = fan(in, pr, eta) 
    
    R = 0.287;   
    
    P01 = in.p0;
    P02 = P01 * pr;
    P1 = in.p;
    P2 = in.p * pr;
    T1 = in.T;
    

    [~,h1,s1,~] = air_props(T1);
    
    fun = @(T2s) entropy_air(T2s) - s1 - R*log(P2/P1); 
    T2s = fzero(fun, T1 * pr^0.3);
    
    [~,h2s,~,~] = air_props(T2s);
    
    h2 = h1 + (h2s - h1) / eta;
    
    fun2 = @(T2) enthalpy_air(T2) - h2;
    T2 = fzero(fun2, T2s);
    
    [~,~,~,gamma2] = air_props(T2);
    h01 = h1 + in.V^2 / 2;
    h02 = h01 + (h2s - h1) / eta;
    out.V = sqrt(2 * (h02 - h2));
    out.M = out.V / sqrt(gamma2*R*T2);

    T02 = T2 * ( 1 + (gamma2 - 1) / 2 * out.M^2 );

    out.T = T2;
    out.p = P2;
    out.T0 = T02;
    out.p0 = P02;
    out.w = h02-h01;
    
    [out.cp,out.h,out.s] = air_props(T2);

end


function out = compressor(in, pr, eta)

  R = 0.287;
  out = in;
  P01 = in.p0;
  P02 = P01 * pr;
  T1 = in.T;
  P1 = in.p;
  P2 = P1*pr;

  [~,h1,s1,~] = air_props(T1);

  fun = @(T2s) entropy_air(T2s) - s1 - R*log(P2/P1);
  T2s = fzero(fun, in.T * (pr^0.3));

  [~,h2s,~,~] = air_props(T2s);

  h2 = h1 + (h2s - h1)/eta;

  fun2 = @(T2) enthalpy_air(T2) - h2;
  T2 = fzero(fun2, T2s);

  [~,~,~,gamma2] = air_props(T2);

  h01 = h1 + in.V^2 / 2;
  h02 = h01 + (h2s - h1) / eta;
  out.V = sqrt(2 * (h02 - h2));
  out.M = out.V / sqrt(gamma2*R*T2);

  T02 = T2 * ( 1 + (gamma2 - 1) / 2 * out.M^2 );

  out.T = T2;
  out.p = P2;
  out.T0 = T02;
  out.p0 = P02;   
  out.dmdt = dmdt_aH;
  out.W = in.mdot*(h02-h01);
  out.p0 = P02;
  out.W = h02-h01;

end

function out = turbine(in, w_req, eta)
  R = 0.287;

  out = in;

  [~,h1,s1,~] = air_props(in.T);

  h01 = h1 + in.V^2 / 2;
  h02 = h01 - eta * w_req;

  V_guess = in.V;
  for i = 1:10
      h2 = h02 - V_guess^2 / 2;
      fun = @(T2) enthalpy_air(T2) - h2;
      out.T = fzero(fun,in.T - 300);

      [~,~,s2,~] = air_props(out.T);
      out.p = in.p * exp((s2 - s1) / R);

      [~,~,~,gamma] = air_props(out.T);
      out.p0 = in.p0 * exp((s2 - s1) / R);
      out.T0 = out.T + V_guess^2 / (2 * (gamma / (gamma - 1)) * R);
      V_new = sqrt(2 * (h02 - h2));

      if abs(V_new - V_guess) < 1e-6
          break
      end
      V_guess = V_new;
  end
out.V = V_guess;
out.M = out.V / sqrt(gamma * R * T2);
out.T0 = out.T + out.V^2 / (2 * (gamma / (gamma - 1)) * R);
out.p0 = in.p0 * exp((s2 - s1) / R); 
out.p  = out.p0 / (1 + (gamma - 1)/2 * out.M^2)^(gamma / (gamma - 1));
[out.cp, out.h, out.s, out.gamma] = air_props(out.T);
end

function out = duct(in, spr)
    out = in;

    out.dmdt = in.dmdt;

    out.p0 = in.p0 * spr;
    out.T0 = in.T0;
end

function out = combustor(in, spr, LHV, eta, R)
    disp(['p0=',num2str(in.p0),', T0=',num2str(in.T0)])
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
    disp(['CEA Input Data: p=',num2str(p*1e-5),' and T=',num2str(T)])
    out.phi = input("φ=");
    out.T_ad = input("T_ad=");
    gamma = input("γ=");
    c = input("c [m/s]=");
    M = 250/c;
    [~, h_reac] = air_props(T);
    out.Q_R = -(input("h_prod=")-h_reac);
    out.p0 = spr*in.p0;
    t0_star = in.T0*(1+gamma*M^2)^2/(2*(gamma+1)*M^2*(1+(gamma-1)/2*M^2));
    out.T0 = t0_star*(in.T0/t0_star + out.Q_R/(gamma*R/(gamma-1)*t0_star));
    out.dmdt_f = out.Q_R/LHV/eta;
end

function out = mixer(core, bypass, spr)
% Assuming 100% mixing takes place
  out = core;
  mdot_core = core.mdot;
  mdot_bypass = bypass.mdot;
  mdot_a = mdot_bypass + mdot_core;

  [~,hc,~,~] = air_props(core.T);
  [~,hb,~,~] = air_props(bypass.T);

  hm = (mdot_core * hc + mdot_bypass * hb) / mdot_a;

  fun = @(T) enthalpy_air(T) - hm;
  T_mix = fzero(fun, (core.T + bypass.T)/2);

  p_mix = min(core.P, bypass.P);
  out.p = p_mix * spr;
  out.T = T_mix;
  out.dmdt = mdot_a;
end

function out = afterburner_inop(in, spr)
    out = duct(in, spr);
end

function out = afterburner_op()
    %TBD probably another combustor
end

function out = nozzle(in, spr, eta, R)
% need to check for choked flow
    out.h0 = in.h0;
    out.p0 = in.p0*spr;
    out.T0 = in.T0;
    i = 0; err = 1; tmp.T = in.T0; tmp.M = 0.5;
    disp(['p0=',num2str(out.p0*1e-5),', T0=',num2str(out.T0)])
    while err > 0.001 && i < 1e4
        i = i+1;
        disp(['i=',num2str(i),', T=',num2str(tmp),', p=',num2str(p)])
        [~, ~, ~, gamma] = air_props(tmp.T);
        out.p = out.p0*(1+(gamma-1)/2*tmp.M^2)^(-gamma/(gamma-1));
        out.u = sqrt(2*eta*gamma/(gamma-1)*R*in.T0*(1-(p/in.p0)^((gamma-1)/gamma)));
        M = out.u/sqrt(gamma*R*tmp.T);
        T = in.T0/(1+(gamma-1)/2*M^2);
        err = sqrt(mean( ((T-tmp.T)/T)^2+((M-tmp.M)/M)^2 ));
        tmp.T = T;
        tmp.M = M;
    end
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
R = 287;
cp_a = air_props(ambient.T);
gamma_a = 1/(1-R/cp_a);
M_a = 0.5;
dmdt_aH = (1-BPR)*dmdt_a;
u = M_a*sqrt(gamma_a*R*ambient.T);
ambient.p0 = ambient.p*(1+(gamma_a-1)/2*M_a^2)^(gamma_a/(gamma_a-1));
ambient.T0 = ambient.T*(1+(gamma_a-1)/2*M_a^2);

engine.diffuser = diffuser(ambient, spr.INT);
engine.fan = fan(engine.diffuser, FAN_pr, eta.FAN);
core = engine.fan;                 %
bypass = engine.fan;               %
core.mdot = dmdt_a / (1 + BPR);    %
bypass.mdot = dmdt_a - core.mdot;  %  Needed for mixer code to run.
bypass = duct(bypass, spr.BPD);    %  
core = duct(core, spr.LPC);        %

engine.lpc = compressor(core, spr.LPC, eta.LPC);
engine.hpc = compressor(engine.lpc, spr.HPC ,eta.HPC);
engine.combustor = combustor(engine.hpc, spr.BRN, LHV, eta.BRN, R);

engine.HPT = turbine(engine.combustor, (engine.lpc.W+engine.fan.W)/core.mdot, eta.HPT);
engine.LPT = turbine(engine.HPT, engine.hpc.W/core.mdot, eta.LPT);
engine.lpc = compressor(engine.fan, spr.LPC, eta.LPC);
engine.hpc = compressor(engine.lpc, spr.HPC ,eta.HPC);
engine.combustor = combustor(engine.hpc, spr.BRN, LHV, eta.BRN, R);
engine.HPT = turbine(engine.combustor, engine.hpc.w/eta.SFT, eta.HPT);
engine.LPT = turbine(engine.HPT, (core.mdot*engine.lpc.w+ambient.dmdt*engine.fan.w)/core.mdot/eta.SFT, eta.LPT);
engine.duct = duct(engine.fan, spr.BPD);
engine.mixer = mixer(engine.LPT, bypass, spr.MXR);
engine.afterburner = afterburner_inop(engine.mixer, spr.ABR);
engine.nozzle = nozzle(engine.afterburner, spr.NOZ, eta.NOZ, R);

% parameters needed
p_4 = engine.hpc.p;
T_4 = engine.hpc.T;
dmdt_f = engine.combustor.dmdt_f;
Q_R = engine.combustor.Q_R;
u_e = engine.nozzle.u;
p_e = engine.nozzle.p;
% derived parameters
f = dmdt_f/dmdt_aH;

% Output Parameters
thrust = (dmdt_a+dmdt_f)*u_e-dmdt_a*u+pi*d_10^2/4*(p_e-p_a); % Thrust, N
ST = thrust/dmdt_a; % Specific thrust, Ns/kg
TSFC = dmdt_f/thrust; % Thrust-specific fuel consumption, kg/Ns

eta_th = (1+f)*(u_e^2-u^2)/(2*f*Q_R); % Thermal efficiency, []
eta_p = thrust*u/dmdt_a/((1+f)*(u_e^2/2)-u^2/2); % Propulsion efficiency, []
eta_0 = eta_th*eta_p; % Overall efficiency, []

LD = 10; % Estimate of cruise lift to drag ratio, []
s = eta_0*LD*Q_R/g*log(1.66); % Aircraft range (assuming fuel is 40% of aircraft weight), m
