function out = analyzePri_oneDevice(topology, f_per, f_sw_typ, I_r, jj, tech, rowID, ganFile, sicFile)
% analyzePri_oneDevice
%
% Analyze primary-side loss and junction temperature for ONE selected switch.
%
% Inputs:
%   topology  : "Multitrack", "2-level Multitrack", or "Fullbridge LLC"
%   f_per     : switching-frequency tolerance, e.g., 0.1
%   f_sw_typ  : typical switching frequency [Hz]
%   I_r       : RMS current [A]
%   jj        : number of parallel switches
%   tech      : "GaN" or "SiC"
%   rowID     : selected row number in the Excel file
%   ganFile   : GaN data file name
%   sicFile   : SiC data file name
%
% Example:
%   out = analyzePri_oneDevice("Multitrack",0.1,500e3,30,4,"GaN",2,...
%                              "GaN Data.xlsx","SiC Data.xlsx");

%% ===================== Defaults =====================
if nargin < 9, sicFile = 'SiC Data.xlsx'; end
if nargin < 8, ganFile = 'GaN Data.xlsx'; end

tech = string(tech);

%% ===================== Column mapping by index =====================
% Same column numbers as your original code

% GaN indices
idx.GaN.Name = 1;
idx.GaN.Rth_jc = 4;
idx.GaN.Cool   = 5;
idx.GaN.Coss   = 7;
idx.GaN.L_min  = 8;
idx.GaN.W_min  = 9;
idx.GaN.Rds    = 10;
idx.GaN.Tjmax  = 11;
idx.GaN.tf     = 13;
idx.GaN.Qg     = 14;
idx.GaN.Vg     = 15;
idx.GaN.Pin_Area = 30;
idx.GaN.R_ds_25 = 26;

% SiC indices
idx.SiC.Name = 1;
idx.SiC.Rth_jc = 4;
idx.SiC.Cool   = 5;
idx.SiC.Coss   = 11;
idx.SiC.L_min  = 12;
idx.SiC.W_min  = 13;
idx.SiC.Rds    = 14;
idx.SiC.Tjmax  = 15;
idx.SiC.tf     = 17;
idx.SiC.Qg     = 18;
idx.SiC.Vg     = 19;
idx.SiC.Pin_Area = 29;

if ~isfield(idx, tech)
    error('tech must be "GaN" or "SiC".');
end

%% ===================== Read selected data table =====================
if tech == "GaN"
    if isempty(ganFile)
        error('ganFile is empty.');
    end
    T = readtable(ganFile,'VariableNamingRule','preserve');
elseif tech == "SiC"
    if isempty(sicFile)
        error('sicFile is empty.');
    end
    T = readtable(sicFile,'VariableNamingRule','preserve');
end

if rowID < 1 || rowID > height(T)
    error('rowID exceeds the number of rows in the selected data file.');
end

map = idx.(tech);
ii = rowID;

%% ===================== Read one selected device =====================
L_min   = T{ii, map.L_min};
W_min   = T{ii, map.W_min};

if isfield(map, 'Pin_Area') && map.Pin_Area <= width(T)
    Pin_Area = T{ii, map.Pin_Area};
    if isempty(Pin_Area) || isnan(Pin_Area)
        Pin_Area = 0;
    end
else
    Pin_Area = 0;
end

Rth_jc   = T{ii, map.Rth_jc};
R_ds_25  = T{ii, map.R_ds_25}; 
fprintf("R_ds_25 is %d\n", R_ds_25);
Coss     = T{ii, map.Coss};
t_f      = T{ii, map.tf};
Q_g      = T{ii, map.Qg};
V_g      = T{ii, map.Vg};
T_j_max  = T{ii, map.Tjmax};
Cooling  = T{ii, map.Cool};

Name = string(T{ii, map.Name});

%% ===================== Known parameters =====================
T_water = 45;

spacing = 4*(0.00254);
copper_oz = 2;
copper_thick = (copper_oz/0.5)*0.00175; % [cm]
board_thick_cm = 0.2;                   % [cm]
Radius_via = 6*(0.00254);

R_via = 0.25*(board_thick_cm)/( ...
    pi*(Radius_via^2) - pi*(Radius_via-copper_thick)^2);

pad_thick_cm = 0.16;

R_plate = 0.08;
Area_plate = 152.4*76.2; % [mm^2]

f_sw_max = f_sw_typ*(1+f_per);

%% ===================== Topology settings =====================
if topology == "Multitrack"
    pri_sw_count = 4;
    duty_cond = 0.5;
elseif topology == "2-level Multitrack"
    pri_sw_count = 8;
    duty_cond = 0.5;
elseif topology == "Fullbridge LLC"
    pri_sw_count = 8;
    duty_cond = 0.5;
else
    error('Unsupported topology in this one-device version.');
end

%% ===================== Current calculation =====================
I_r_pk = I_r*sqrt(2);
I_rating = I_r_pk;

% Current through each parallel device
I_d = I_rating / jj;

%% ===================== Area calculation =====================
Area_single_device = L_min * W_min + Pin_Area;
Area_total = pri_sw_count * jj * Area_single_device;

%% ===================== Thermal resistance calculation =====================
N_L = floor(L_min*0.1/(2*Radius_via + 2*spacing));
N_W = floor(W_min*0.1/(2*Radius_via + 2*spacing));
N_vias_max = N_L*N_W;

Area_vias = N_vias_max*pi*((Radius_via*10)^2); % [mm^2]
Area_fr4  = L_min*W_min - Area_vias;           % [mm^2]

if Area_fr4 <= 0
    warning('Area_fr4 is non-positive. Check via and footprint dimensions.');
    Area_fr4 = eps;
end

R_fr4 = 4350*(board_thick_cm*10)/Area_fr4;

if string(Cooling) == "Top"
    Rth_board_min = 0;
else
    if N_vias_max <= 0
        warning('No vias fit in the selected footprint. Setting Rth_board_min = R_fr4.');
        Rth_board_min = R_fr4;
    else
        Rth_board_min = ((R_via/N_vias_max)*R_fr4) / ...
                        ((R_via/N_vias_max) + R_fr4);
    end
end

Rth_pw    = R_plate * Area_plate / (L_min * W_min);
Rth_inter = 6*pad_thick_cm / (L_min * W_min * 0.01);

Rth_total_jw = Rth_jc + Rth_inter + Rth_board_min + Rth_pw;

%% ===================== Rds,on temperature equation =====================
% Normalized Rds,on curve:
% Tj in degC
% Rds_norm is normalized to the datasheet curve
Rds_norm = @(Tj) -8.3995e-10*Tj.^4 ...
                 +1.3410e-7*Tj.^3 ...
                 +2.1477e-5*Tj.^2 ...
                 +4.7680e-3*Tj ...
                 +0.86269;

% Convert normalized curve to actual Rds,on.
% R_ds_25 is the table value at 25 degC.
Rds_on = @(Tj) R_ds_25 .* Rds_norm(Tj)./ Rds_norm(25);

%% ===================== Iterative loss-temperature calculation =====================
max_iter = 100;
tol = 1e-6;

T_j_old = T_water + 20; % initial guess

for iter = 1:max_iter

    % Step 1: calculate Rds,on based on current temperature guess
    R_ds_now = Rds_on(T_j_old);

    % Step 2: calculate loss using this Rds,on
    P_cond_position = duty_cond * jj * I_d^2 * R_ds_now;

    P_off_position = jj * ((I_d^2) * ((t_f*1e-9)^2) / ...
        (24*2*Coss*1e-12)) / (1/f_sw_max);

    P_gate_position = jj * V_g * Q_g * 1e-9 * f_sw_max;

    % Loss of one switch position
    P_position = P_cond_position + P_off_position + P_gate_position;

    % Loss per physical device
    P_per_device = P_position / jj;

    % Step 3: update junction temperature based on loss
    T_j_new = T_water + P_per_device * Rth_total_jw;

    % Step 4: check convergence
    if abs(T_j_new - T_j_old) < tol
        break;
    end

    T_j_old = T_j_new;
end

T_j = T_j_new;
R_ds_final = Rds_on(T_j);

%% ===================== Final total losses =====================
P_cond_total = pri_sw_count * P_cond_position;
P_off_total  = pri_sw_count * P_off_position;
P_gate_total = pri_sw_count * P_gate_position;
P_total      = P_cond_total + P_off_total + P_gate_total;

thermal_pass = T_j < (T_j_max - 30);

%% ===================== Output =====================
out = struct();

out.Tech = tech;
out.Name = Name;
out.RowID = rowID;
out.Topology = topology;
out.jj = jj;

out.Area_single_device = Area_single_device;
out.Area_total = Area_total;

out.Rds_25 = R_ds_25;
out.Rds_final = R_ds_final;
out.T_j = T_j;
out.T_j_max = T_j_max;
out.ThermalPass = thermal_pass;
out.Iterations = iter;

out.P_cond_total = P_cond_total;
out.P_off_total = P_off_total;
out.P_gate_total = P_gate_total;
out.P_total = P_total;
out.P_per_device = P_per_device;

out.Rth_jc = Rth_jc;
out.Rth_board = Rth_board_min;
out.Rth_interface = Rth_inter;
out.Rth_plate = Rth_pw;
out.Rth_total_jw = Rth_total_jw;

out.N_vias = N_vias_max;
out.R_via_single = R_via;
out.R_fr4 = R_fr4;

%% ===================== Display result =====================
fprintf('\n========== One Device Result ==========\n');
fprintf('Device: %s\n', Name);
fprintf('Technology: %s\n', tech);
fprintf('Row ID: %d\n', rowID);
fprintf('Parallel number jj: %d\n', jj);
fprintf('Final Rds,on: %.6g Ohm\n', R_ds_final);
fprintf('Final Tj: %.2f degC\n', T_j);
fprintf('Total loss: %.3f W\n', P_total);
fprintf('  Conduction loss: %.3f W\n', P_cond_total);
fprintf('  Turn-off loss: %.3f W\n', P_off_total);
fprintf('  Gate loss: %.3f W\n', P_gate_total);
fprintf('Iterations: %d\n', iter);


fprintf('\n--- New code thermal check ---\n');
fprintf('R_ds_25 = %.9f Ohm\n', R_ds_25);
fprintf('R_ds_final = %.9f Ohm\n', R_ds_final);
fprintf('P_per_device = %.6f W\n', P_per_device);
fprintf('Rth_jc = %.6f\n', Rth_jc);
fprintf('Rth_inter = %.6f\n', Rth_inter);
fprintf('Rth_board_min = %.6f\n', Rth_board_min);
fprintf('Rth_pw = %.6f\n', Rth_pw);
fprintf('Rth_total_jw = %.6f\n', Rth_total_jw);
fprintf('T_water = %.6f\n', T_water);
fprintf('T_j = %.6f\n', T_j);
fprintf('T_j_max = %.6f\n', T_j_max);
fprintf('Limit = %.6f\n', T_j_max - 30);
fprintf('P_cond = %.6f\n', P_cond_position);
fprintf('P_off = %.6f\n', P_off_position);
fprintf('P_gate = %.6f\n', P_gate_position);


if thermal_pass
    fprintf('Thermal check: PASS\n');
else
    fprintf('Thermal check: FAIL\n');
end

end