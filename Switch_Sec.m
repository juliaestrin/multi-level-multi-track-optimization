clear
SecData = readtable('Sec Data.xlsx');

%% Known Parameters
Vin = 1500;
Vin_per = 0.1;
Vin_max = Vin*(1+Vin_per);
Vout = 48; 
R_via = 134.83; % [C/W]
spacing = 4*(0.00254); % [cm]
D_via = 6*(0.00254); % [cm]
T_water = 45;
R_plate = 0.08;
Area_plate = 152.4*76.2; % [mm2]
%R_interface = 1.01;
f_sw_typ= 1000e3;
f_sw_max = f_sw_typ*(1.25);
f_sw_min = f_sw_typ*(0.75);
Power = 6.25e3;
N_tr = 2;
N_x = 4;
N_s = 1;
R_load = (Vout^2)/Power;

%% Calculate the I_rms
I_sec_pk = (pi/2)*(Vout/R_load); % I_oe*8*sqrt(2)
I_rating = ceil(I_sec_pk/5)*5; % Add some offset? 

%% Calculate the Loss for each switch
% ROW: different switch; COLUMN: number of parallel switch
data_size = size(SecData);
n_sw = data_size(1);
%n_sw=1;
Area = nan(n_sw,10);
P_sw1 = nan(n_sw,10);
P_sw2 = nan(n_sw,10);
P_sw3 = nan(n_sw,10);
P_cond = nan(n_sw,10);
P_total = nan(n_sw,10);
P_gate = nan(n_sw, 10);
P_total_plot = nan(n_sw,10);
T_c = nan(n_sw,10);
T_j = nan(n_sw,10);
T_j_plot = nan(n_sw,10);
T_per = nan(n_sw,10);
Loss_Area = nan(n_sw,10);

for ii=1:n_sw
    L_min = SecData{ii,7}; % [mm]
    W_min = SecData{ii,8}; % [mm]
    N_L = floor(L_min*0.1/(D_via+2*spacing));
    N_W = floor(W_min*0.1/(D_via+2*spacing));
    N_vias_max = N_L*N_W;
    Area_vias = N_vias_max*pi*((D_via*10/2)^2); % [mm2]
    Area_fr4 = L_min*W_min-Area_vias; % [mm2]
    R_fr4 = 4350*1.6/Area_fr4;
    %R_fr4 = 69.9;
    %fprintf("The N_vias_max is %d\n", N_vias_max);
    %fprintf("The R_fr4 is %d\n", R_fr4);
    Rth_board_min = ((R_via/N_vias_max)*R_fr4/((R_via/N_vias_max)+R_fr4));
    Rth_jc = SecData{ii,4};
    R_ds_max = SecData{ii,9};
    Coss = SecData{ii,5}; % [pF]
    t_r = SecData{ii,11}; % [ns]
    t_f = SecData{ii,12}; % [ns]
    Q_g = SecData{ii,13}; % [nC]
    V_g = SecData{ii,14};

    for jj=1:10
        Area(ii,jj) = (jj)*L_min*W_min; % [mm2]P
        I_d = I_rating/(jj); % the pk current through each switch
        P_cond(ii,jj) = 0.5*(jj)*(I_d^2)*R_ds_max;
        P_gate(ii,jj) = (jj)*V_g*Q_g*(1e-9)*f_sw_max;
        P_total(ii,jj) = P_cond(ii,jj) + 0*P_gate(ii,jj);
        % Plate-to-water thermal resistance
        Rth_pw = R_plate*Area_plate/(1*L_min*W_min);
        % Interface thermal resistance
        Rth_inter = 6*(0.6)/(L_min*W_min*0.01);
        % fprintf("The Rth_inter is %d\n", Rth_inter);
        T_c(ii,jj) = (P_total(ii,jj)/(jj))*(Rth_inter + Rth_board_min + Rth_pw) + T_water;
        T_j(ii,jj) = (P_total(ii,jj)/(jj))*(Rth_jc + Rth_inter + Rth_board_min + Rth_pw) + T_water;
        T_j_max = SecData{ii,10};
        Loss_Area(ii,jj) = Area(ii,jj)*P_total(ii,jj);

        if T_j(ii,jj) < (T_j_max - 30)
            T_per(ii,jj) = (T_j_max - T_j(ii,jj))/T_j_max;
            P_total_plot(ii,jj) = P_total(ii,jj);
            T_j_plot(ii,jj) = T_j(ii,jj);
        else
            P_total_plot(ii,jj) = NaN;
            T_j_plot(ii,jj) = NaN;
        end
    end
end

%% Setting the Colors
colors = [
    0.1216    0.4667    0.7059   % blue
    1.0000    0.4980    0.0549   % orange
    0.1725    0.6275    0.1725   % green
    0.8392    0.1529    0.1569   % red
    0.5804    0.4039    0.7412   % purple
    0.5490    0.3373    0.2941   % brown
    0.8902    0.4667    0.7608   % pink
    0.4980    0.4980    0.4980   % gray
    0.7373    0.7412    0.1333   % yellow-green
    0.0902    0.7451    0.8118   % teal
    0.7373    0.2235    0.2235   % dark red
    0.2549    0.4118    0.8824   % bright blue
    0.4940  0.1840  0.5560   % dark purple
    0.9290  0.6940  0.1250   % gold
    0.3010  0.7450  0.9330   % cyan
];

% colors = [
%     0.1216  0.4667  0.7059   % blue
%     1.0000  0.4980  0.0549   % orange
%     0.1725  0.6275  0.1725   % green
%     0.8392  0.1529  0.1569   % red
%     0.5804  0.4039  0.7412   % purple
%     0.5490  0.3373  0.2941   % brown
%     0.8902  0.4667  0.7608   % pink
%     0.4980  0.4980  0.4980   % gray
%     0.7373  0.7412  0.1333   % olive
%     0.0902  0.7451  0.8118   % teal
%     0.7373  0.2235  0.2235   % dark red
%     0.2549  0.4118  0.8824   % bright blue
%     0.0000  0.0000  0.0000   % black
%     0.0000  0.6000  0.5000   % blue-green
%     0.6000  0.6000  0.0000   % dark yellow
%     0.0000  0.4470  0.7410   % navy
%     0.8500  0.3250  0.0980   % dark orange
%     0.9290  0.6940  0.1250   % gold
%     0.4940  0.1840  0.5560   % dark purple
%     0.3010  0.7450  0.9330   % cyan
% ];

%% Plot the Graphs
figure(1)
% Plot the total power loss
t = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile;
hold on;
for i = 1:n_sw
    plot(1:10, P_total_plot(i, :),'-o','Color', ...
        colors(i,:),'MarkerFaceColor', colors(i,:), ...
        'LineWidth', 1.5,'MarkerSize', 6);   % Plot each row
end
xlabel('Number of Parallel Switches','interpreter','latex','FontSize',15);
ylabel('Total Power Loss [W]','interpreter','latex','FontSize',15);
%legend(SecData.Name, 'Location', 'best','FontSize',10);
grid on;

% Plot the T_per
nexttile;
hold on; 
for i = 1:n_sw
    plot(1:10, T_j_plot(i, :),'-o','Color', ...
        colors(i,:),'MarkerFaceColor', colors(i,:), ...
        'LineWidth', 1.5,'MarkerSize', 6);   % Plot each row
end
xlabel('Number of Parallel Switches','interpreter','latex','FontSize',15);
ylabel('Junction Temperature [$^{\circ}\mathrm{C}$]','interpreter','latex','FontSize',15);
legend(SecData.Name, 'Location', 'best','FontSize',10);
sgtitle('Si/GaN MOSFET with $V_{ds} = 100$ V and $I_d \geq 90$ A','interpreter','latex','FontSize',15);
grid on;













