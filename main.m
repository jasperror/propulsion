% ENGR477 Propulsion
% Project 1: F135 Engine
% Jasper Palmer & Jackson Gilbert

clear; close all; clc;

%% Constants

% Design-point Operating Conditions
dmdt_tot = 150; % Total air mass flow rate, kg/s
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


%% Inlet ambiant conditions 
gamma_guess = 1.4;
flow.P0 = p_a * (1 + ((gamma_guess - 1)/2) * M_1^2)^(gamma_guess/(gamma_guess - 1)); % Ambiant static pressure, Pa
flow.T0 = T_a * (1 + ((gamma_guess - 1)/2) * M_1^2); % Ambiant static temperature, K
flow.mdot = dmdt_tot; % Mass flow rate, kg/s
flow.M = M_1;
[flow.cp, flow.h, flow.s] = air_props(flow.T0);
%% Engine Components

function [cp,h,s] = air_props(T)

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

end

function s = entropy_air(T)
    [~,~,s] = air_props(T);
end

function h = enthalpy_air(T)
    [~,h,~] = air_props(T);
end
function out = diffuser(in,pr) %Isentropic and adiabatic
out = in;
out.P0 = in.P0 * pr; 

% adiabadic -> T0_a = T0_2

end

function out = fan(in, pr, eta) 

R = 0.287;   

P1 = in.P0;
P2 = P1 * pr;

[~,h1,s1] = air_props(in.T0);

fun = @(T2s) entropy_air(T2s) - s1 - R*log(P2/P1);
T2s = fzero(fun, in.T0 * pr^0.3);

[~,h2s,~] = air_props(T2s);

h2 = h1 + (h2s - h1)/eta;

fun2 = @(T2) enthalpy_air(T2) - h2;
T2 = fzero(fun2, T2s);

out.T0 = T2;
out.P0 = P2;

[out.cp,out.h,out.s] = air_props(T2);

end

function out = compressor(in, pr, eta)

R = 0.287;
out = in;
P1 = in.P0;
P2 = P1 * pr;


[~,h1,s1] = air_props(in.T0);


fun = @(T2s) entropy_air(T2s) - s1 - R*log(P2/P1);
T2s = fzero(fun, in.T0 * (pr^0.3));

[~,h2s,~] = air_props(T2s);

h2 = h1 + (h2s - h1)/eta;

fun2 = @(T2) enthalpy_air(T2) - h2;
T2 = fzero(fun2, T2s);

out.T0 = T2;
out.P0 = P2;  
[out.cp,out.h,out.s] = air_props(T2);

end


function out = combustor()

end

function out = turbine(in, W_req, eta)

R = 0.287;

out = in;

[~,h1,s1] = air_props(in.T0);

h2 = h1 - W_req;

fun = @(T2) enthalpy_air(T2) - h2;
T2 = fzero(fun, in.T0 - 300);

[~,~,s2] = air_props(T2);

P2 = in.P0 * exp((s2 - s1)/R);

P2 = in.P0 * (P2/in.P0)^eta;

out.T0 = T2;
out.P0 = P2;

[out.cp,out.h,out.s] = air_props(T2);

end

function out = duct(in, spr)
out = in;
out.P0 = in.P0 * spr;
out.T0 = in.T0;
[out.cp,out.h,out.s] = air_props(out.T0);

end

function out = mixer()
% Assuming 100% mixing takes place

end

function out = afterburner()

end

function out = nozzle()
% Nozzle type C
end
%% Engine Structure

%{
Element Numbering
a: Atmospheric (inlet condition)
1: Fan inlet
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