% Draw 3D Figure of EQ Core Transformer Geometry
%
% Author:  Julia Estrin
% Date:    02-12-2026
% Updated: 03-27-2026 - Added support for round vs square centerpost
%
% Description:
%   Renders a 3D visualization of the EQ-core transformer geometry using
%   dimensions from the result struct produced by calculate_VIRT_design.
%   Core sections (yokes, legs, center post) and primary winding are drawn
%   as colored solid patches. All dimensions are displayed in mm.
%   Supports both round (cylindrical) and square centerpost geometries.
%
% Inputs:
%   result    - Design result struct from calculate_VIRT_design
%               Must include 'centerpost_shape' field ('round' or 'square')
%   Pv_max    - Maximum core loss density [W/m^3]  (used for figure title)
%   w_height  - Winding window height [m]           (used for figure title)
%   material  - Material name string (used for figure title)
%   T_tx      - Transformer temperature [°C] (used for figure title)
%
% Usage:
%   core3Dfigure(result, Pv_max, w_height, material, T_tx);

function core3Dfigure(result, Pv_max, w_height, material, T_tx)

    figure(1);
    clf;
    hold on;
    axis equal;
    view(35, 15);
    grid on;
    grid minor;
    box on;

    % --- Colors ---
    core_color = [0.5, 0.5, 0.8];
    copper_pri = [0.85, 0.45, 0.01];

    % --- Determine centerpost shape ---
    if isfield(result, 'centerpost_shape')
        centerpost_shape = result.centerpost_shape;
    else
        % Default to square for backward compatibility
        centerpost_shape = 'square';
    end

    % --- Convert all dimensions to mm ---
    l_core    = result.l_core    * 1e3;
    w_core    = result.w_core    * 1e3;
    l_leg     = result.l_leg     * 1e3;
    h_yoke    = result.h_yoke    * 1e3;
    w_b       = result.w_b       * 1e3;
    h_core    = result.h_core    * 1e3;
    l_winding = result.l_winding * 1e3;
    w_winding = result.w_winding * 1e3;
    t_pcb     = 1.6;                        % [mm] PCB thickness (fixed)

    % Get centerpost radius (for round) or compute from w_core (for square)
    if strcmp(centerpost_shape, 'round') && isfield(result, 'r_centerpost') && ~isnan(result.r_centerpost)
        r_centerpost = result.r_centerpost * 1e3;
    else
        r_centerpost = w_core / 2;
    end

    % --- Bottom Yoke ---
    drawBox([-l_core/2, 0, 0], [l_core, w_core, h_yoke], core_color);

    % --- Top Yoke ---
    z_top = h_yoke + w_b;
    drawBox([-l_core/2, 0, z_top], [l_core, w_core, h_yoke], core_color);

    % --- Left Leg ---
    drawBox([-l_core/2, 0, h_yoke], [l_leg, w_core, w_b], core_color);

    % --- Right Leg ---
    x_right = l_core - l_leg;
    drawBox([x_right - l_core/2, 0, h_yoke], [l_leg, w_core, w_b], core_color);

    % --- Center Post (shape-dependent) ---
    if strcmp(centerpost_shape, 'round')
        % Cylindrical center post
        % Spans from bottom yoke top to top yoke bottom (height = w_b)
        x_center = l_core / 2;
        [X, Y, Z] = cylinder(r_centerpost, 50);  % 50 facets for smooth cylinder
        Z = Z * w_b;
        X = X + x_center - l_core/2;
        Y = Y + w_core / 2;
        surf(X, Y, Z + h_yoke, 'FaceColor', core_color, 'EdgeColor', 'none');
        
        % Add top and bottom caps for the cylinder
        theta = linspace(0, 2*pi, 51);
        x_cap = r_centerpost * cos(theta);
        y_cap = r_centerpost * sin(theta) + w_core/2;
        
        % Bottom cap
        fill3(x_cap, y_cap, ones(size(x_cap)) * h_yoke, core_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
        % Top cap
        fill3(x_cap, y_cap, ones(size(x_cap)) * (h_yoke + w_b), core_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
    else
        % Square center post (box)
        drawBox([-w_core/2, 0, h_yoke], [w_core, w_core, w_b], core_color);
    end

    % --- Primary Winding (4 PCB strips forming a rectangular frame) ---
    % Winding sits at mid-height of the window, centered at z = h_yoke + w_b/2
    z_winding = h_yoke + w_b/2 - t_pcb/2;

    % Front strip (negative y side)
    drawBox([l_leg - l_core/2,  -(w_winding/2 - w_core/2),  z_winding], ...
            [l_winding,           w_winding/2 - r_centerpost, t_pcb], copper_pri);

    % Back strip (positive y side)
    drawBox([l_leg - l_core/2,   w_core,                     z_winding], ...
            [l_winding,           w_winding/2 - r_centerpost, t_pcb], copper_pri);

    % Left strip
    drawBox([l_leg - l_core/2,  -(w_winding/2 - w_core/2),  z_winding], ...
            [l_winding/2 - r_centerpost, w_winding,           t_pcb], copper_pri);

    % Right strip
    drawBox([l_core/2 + r_centerpost - l_core/2, -(w_winding/2 - w_core/2), z_winding], ...
            [l_winding/2 - r_centerpost,           w_winding,                t_pcb], copper_pri);

    % --- Formatting ---
    xlabel('X [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Z [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Build title string - include centerpost shape info
    title_str = sprintf(['3D EQ-Core Transformer (%s centerpost) - %s\n' ...
                         'P_v = %.0f kW/m³,  Window Height = %.1f mm (x-axis),  Total Loss = %.2f W\n' ...
                         'Copper Loss = %.2f W,  Core Loss = %.2f W\n' ...
                         'Temp = %.2f °C'], ...
                        centerpost_shape, ...
                        material, ...
                        Pv_max / 1e3, ...
                        w_height * 1e3, ...
                        result.P_total, ...
                        result.P_pri + result.P_sec, ...
                        result.P_core, ...
                        T_tx);
    title(title_str, 'FontSize', 13, 'FontWeight', 'bold');

    rotate3d on;
    hold off;

end

% % Draw 3D Figure of EQ Core Transformer Geometry
% %
% % Author:  Julia Estrin
% % Date:    02-12-2026
% %
% % Description:
% %   Renders a 3D visualization of the EQ-core transformer geometry using
% %   dimensions from the result struct produced by calculate_VIRT_design.
% %   Core sections (yokes, legs, center post) and primary winding are drawn
% %   as colored solid patches. All dimensions are displayed in mm.
% %
% % Inputs:
% %   result    - Design result struct from calculate_VIRT_design
% %   Pv_max    - Maximum core loss density [W/m^3]  (used for figure title)
% %   w_height  - Winding window height [m]           (used for figure title)
% %
% % Usage:
% %   core3Dfigure(result, Pv_max, w_height);
% 
% function core3Dfigure(result, Pv_max, w_height, material, T_tx)
% 
%     figure(1);
%     clf;
%     hold on;
%     axis equal;
%     view(35, 15);
%     grid on;
%     grid minor;
%     box on;
% 
%     % --- Colors ---
%     core_color = [0.5, 0.5, 0.8];
%     copper_pri = [0.85, 0.45, 0.01];
% 
%     % --- Convert all dimensions to mm ---
%     l_core        = result.l_core        * 1e3;
%     w_core        = result.w_core        * 1e3;
%     l_leg         = result.l_leg         * 1e3;
%     h_yoke        = result.h_yoke        * 1e3;
%     w_b           = result.w_b           * 1e3;
%     h_core        = result.h_core        * 1e3;
%     %r_centerpost  = result.r_centerpost  * 1e3;
%     l_winding     = result.l_winding     * 1e3;
%     w_winding     = result.w_winding     * 1e3;
%     t_pcb         = 1.6;                            % [mm] PCB thickness (fixed)
%     r_centerpost = w_core/2;
% 
%     % --- Bottom Yoke ---
%     drawBox([-l_core/2, 0, 0], [l_core, w_core, h_yoke], core_color);
% 
%     % --- Top Yoke ---
%     z_top = h_yoke + w_b;
%     drawBox([-l_core/2, 0, z_top], [l_core, w_core, h_yoke], core_color);
% 
%     % --- Left Leg ---
%     drawBox([-l_core/2, 0, h_yoke], [l_leg, w_core, w_b], core_color);
% 
%     % --- Right Leg ---
%     x_right = l_core - l_leg;
%     drawBox([x_right - l_core/2, 0, h_yoke], [l_leg, w_core, w_b], core_color);
% 
%     % --- Center Post (cylinder along z) ---
%     % Spans from bottom yoke top to top yoke bottom (height = w_b)
%     % x_center = l_core / 2;
%     % [X, Y, Z] = cylinder(r_centerpost);
%     % Z = Z * w_b;
%     % X = X + x_center - l_core/2;
%     % Y = Y + w_core / 2;
%     % surf(X, Y, Z + h_yoke, 'FaceColor', core_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
% 
%     % Center Leg 
%     drawBox([-w_core/2, 0, h_yoke], [w_core, w_core, w_b], core_color);
% 
%     % --- Primary Winding (4 PCB strips forming a rectangular frame) ---
%     % Winding sits at mid-height of the window, centered at z = h_yoke + w_b/2
%     z_winding = h_yoke + w_b/2 - t_pcb/2;
% 
%     % Front strip (negative y side)
%     drawBox([l_leg - l_core/2,  -(w_winding/2 - w_core/2),  z_winding], ...
%             [l_winding,           w_winding/2 - r_centerpost, t_pcb], copper_pri);
% 
%     % Back strip (positive y side)
%     drawBox([l_leg - l_core/2,   w_core,                     z_winding], ...
%             [l_winding,           w_winding/2 - r_centerpost, t_pcb], copper_pri);
% 
%     % Left strip
%     drawBox([l_leg - l_core/2,  -(w_winding/2 - w_core/2),  z_winding], ...
%             [l_winding/2 - r_centerpost, w_winding,           t_pcb], copper_pri);
% 
%     % Right strip
%     drawBox([l_core/2 + r_centerpost - l_core/2, -(w_winding/2 - w_core/2), z_winding], ...
%             [l_winding/2 - r_centerpost,           w_winding,                t_pcb], copper_pri);
% 
%     % --- Formatting ---
%     xlabel('X [mm]', 'FontSize', 12, 'FontWeight', 'bold');
%     ylabel('Y [mm]', 'FontSize', 12, 'FontWeight', 'bold');
%     zlabel('Z [mm]', 'FontSize', 12, 'FontWeight', 'bold');
%     title(sprintf(['3D EE-Core Transformer - ', material, '\n' ...
%                    'P_v = %.0f kW/m³,  Window Height = %.1f mm (x-axis),  Total Loss = %.2f W\n' ...
%                    'Copper Loss = %.2f W,  Core Loss = %.2f W \n Temp = = %.2f °C'], ...
%             Pv_max / 1e3, w_height * 1e3, ...
%            result.P_total, result.P_pri + result.P_sec, result.P_core, T_tx), ...
%            'FontSize', 13, 'FontWeight', 'bold');
% 
%     rotate3d on;
%     hold off;
% 
% end