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

%% Engine Components

function out = diffuser()

end

function out = nozzle()
% Nozzle type C
end

function out = fan()

end

function out = compressor()

end

function out = turbine()

end

function out = combustor()

end

function out = duct()

end

function out = mixer()
% Assuming 100% mixing takes place

end

function out = afterburner()

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