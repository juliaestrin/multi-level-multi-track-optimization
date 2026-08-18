% Check specific device and number of parallel device

clear 
addpath('Switch Functions');

%% Device
% GaN GS065-LR
Data = readtable("Switch Functions\OneDevice.xlsx");

spacing = 4*(0.00254);
via_pad = 5*(0.00254); % 4-6 mil 
copper_oz = 2;
copper_thick = (copper_oz/0.5)*0.00175; % [cm]
board_thick_cm = 0.2472; % [cm]
Radius_via = 6*(0.00254); % [cm]
R_via = 0.25*(board_thick_cm)/(pi*(Radius_via^2)-pi*(Radius_via-copper_thick)^2);
%R_via = 0.25*(board_thick_cm)/(pi*(Radius_via^2));
fprintf('the R_via is %d\n', R_via);

% Interface thickness 
inter_thick_cm = 0.16; % [cm]

% Cooling water temperature
T_water = 45;

% Frequency range
f_sw_typ = 1000e3;
f_per = 0.25;
f_sw_max = f_sw_typ*(1+f_per);

% Number of parallel device
jj = 1;

% Device Parameter
Rth_jc   = Data{1,4};
L_min    = Data{1,8};
W_min    = Data{1,9};
R_ds_max = Data{1,10};
Coss     = Data{1,7};
t_f      = Data{1,13};
Q_g      = Data{1,14};
V_g      = Data{1,15};

%% Thermal Resistance Calculation

% Thermal Vias
Thermal_L = Data{1,26};
Thermal_W = Data{1,27};
N_L = floor(Thermal_L*0.1/(2*Radius_via+2*spacing+2*via_pad));
N_W = floor(Thermal_W*0.1/(2*Radius_via+2*spacing+2*via_pad));

N_vias_max = N_L*N_W;
fprintf("the maximum number of vias is %d \n", N_vias_max);

Area_vias = N_vias_max*pi*((Radius_via*10)^2); % [mm2]
Area_fr4  = Thermal_L*Thermal_W - Area_vias; % [mm2]
R_fr4     = 4350*(board_thick_cm*10)/Area_fr4;
fprintf("the R_fr4 is %d\n", R_fr4);

Rth_board_min = ((R_via/N_vias_max)*R_fr4/((R_via/N_vias_max)+R_fr4));
fprintf("the R_board is %d\n",Rth_board_min);

Rth_pw    = 1; % Fix this to be one

Rth_inter = (1/0.178)*(inter_thick_cm) / (Thermal_L* Thermal_W * 0.01); % [cm]
fprintf("the R_inter is %d\n", Rth_inter);

%% Loss Calculation
Ir_pk = 13.8563;
I_d = Ir_pk/jj;

% Conduction Loss
duty_cond = 0.5;
P_cond = duty_cond * jj * (I_d^2) * R_ds_max;

% Turn-off Loss
P_off  = jj * ((I_d^2) * ((t_f*1e-9)^2) / (24*2*Coss*1e-12)) / (1/f_sw_max);

% Gate Loss
P_gate = jj * V_g * Q_g * 1e-9 * f_sw_max;

P_total = P_off + P_cond + P_gate;

T_j = (P_total/jj) * (Rth_jc   + Rth_inter + Rth_board_min + Rth_pw) + T_water;

fprintf("the final T_j is %d\n", T_j);
