% Create Design Summary Figure
%
% Author:  Julia Estrin
% Date:    03-27-2026
%
% Description:
%   Generates a professional summary figure displaying all key design
%   specifications for the Multilevel Multitrack converter.
%   Box heights are automatically calculated based on content.
%
% Usage:
%   design_summary_figure(topology, Vin, Vo, Pmax, fsw, f0, nt, ...
%       LLC_design, opt, material, stackup, t_cu_pri, t_cu_sec, T_tx);

function fig = design_summary_figure(topology, Vin, Vo, Pmax, fsw, f0, nt, ...
    LLC_design, opt, material, stackup, t_cu_pri, t_cu_sec, T_tx)

    % Create figure
    fig = figure();

    % Colors
    headerColor = [0.2, 0.3, 0.5];
    valueColor = [0.1, 0.1, 0.1];
    labelColor = [0.35, 0.35, 0.35];
    boxColor = [0.95, 0.95, 0.98];
    accentColor = [0.85, 0.33, 0.1];

    % Create main axes
    ax = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
    hold on;

    % =====================================================================
    % LAYOUT PARAMETERS
    % =====================================================================
    dy = 0.036;             % Line spacing (vertical distance between lines)
    header_pad = 0.045;     % Padding from top of box to first line
    bottom_pad = 0.025;     % Padding from last line to bottom of box
    
    % Calculate box height based on number of lines
    calc_height = @(n_lines) header_pad + (n_lines * dy) + bottom_pad;
    
    % Column positions and widths
    col_gap = 0.02;
    col_w = 0.30;
    col1_x = 0.025;
    col2_x = col1_x + col_w + col_gap;
    col3_x = col2_x + col_w + col_gap;
    
    % Row 1 starts below title
    row1_top = 0.89;
    
    % Content line counts for each section
    n_electrical = 7;
    n_llc = 7;
    n_material = 4;  % includes material name
    n_winding = 2;
    n_geometry = 7;  % includes centerpost label
    n_losses = 5;    % includes highlighted total
    n_thermal = 3;
    
    % Calculate section heights
    h_electrical = calc_height(n_electrical);
    h_llc = calc_height(n_llc);
    h_material = calc_height(n_material);
    h_winding = calc_height(n_winding);
    h_geometry = calc_height(n_geometry);
    h_losses = calc_height(n_losses) + 0.02;  % extra for highlight box
    h_thermal = calc_height(n_thermal);
    
    % Row 2 top position (below tallest box in row 1 + gap)
    row1_max_h = max([h_electrical, h_llc, h_material + h_winding + col_gap]);
    row_gap = 0.03;
    row2_top = row1_top - row1_max_h - row_gap;

    % =====================================================================
    % TITLE
    % =====================================================================
    text(0.5, 0.965, 'TRANSFORMER/LLC DESIGN SUMMARY', ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', headerColor, ...
        'HorizontalAlignment', 'center');
    
    text(0.5, 0.93, sprintf('%s Topology', topology), ...
        'FontSize', 11, 'Color', labelColor, ...
        'HorizontalAlignment', 'center');

    % =====================================================================
    % SECTION 1: ELECTRICAL (Top Left)
    % =====================================================================
    sx = col1_x; sy = row1_top; sw = col_w; sh = h_electrical;
    draw_section_box(sx, sy - sh, sw, sh, boxColor);
    draw_header(sx + sw/2, sy - 0.022, 'ELECTRICAL', headerColor);
    
    y = sy - header_pad;
    xl = sx + 0.012; xv = sx + sw - 0.012;
    
    draw_line(xl, xv, y, 'Input Voltage', sprintf('%.0f V', Vin), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'Output Voltage', sprintf('%.0f V', Vo), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'Output Power', sprintf('%.2f kW', Pmax/1e3), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'Tracks (series)', sprintf('%d', nt), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'Turns Ratio', sprintf('%d : 1/2', LLC_design.N), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'Switching Freq', sprintf('%.0f kHz', fsw/1e3), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'Core Freq', sprintf('%.2f MHz', f0/1e6), labelColor, valueColor);

    % =====================================================================
    % SECTION 2: LLC RESONANT TANK (Top Center)
    % =====================================================================
    sx = col2_x; sy = row1_top; sw = col_w; sh = h_llc;
    draw_section_box(sx, sy - sh, sw, sh, boxColor);
    draw_header(sx + sw/2, sy - 0.022, 'LLC RESONANT TANK', headerColor);
    
    y = sy - header_pad;
    xl = sx + 0.012; xv = sx + sw - 0.012;
    
    draw_line_sub(xl, xv, y, 'L', 'r', '', sprintf('%.3f μH', LLC_design.Lr*1e6), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'L', 'm', '', sprintf('%.3f μH', LLC_design.Lm*1e6), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'C', 'r', '', sprintf('%.3f μF', LLC_design.Cr*1e6), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'I', 'r,rms', '', sprintf('%.2f A', LLC_design.Ir_rms), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'I', 'r,pk', '', sprintf('%.2f A', LLC_design.Ir_rms*sqrt(2)), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'Q', 'e,max', '', sprintf('%.4f', LLC_design.Qe_max), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'L', 'n', '', sprintf('%.1f', LLC_design.Lm/LLC_design.Lr), labelColor, valueColor);

    % =====================================================================
    % SECTION 3: CORE MATERIAL (Top Right - Upper)
    % =====================================================================
    sx = col3_x; sy = row1_top; sw = col_w; sh = h_material;
    draw_section_box(sx, sy - sh, sw, sh, boxColor);
    draw_header(sx + sw/2, sy - 0.022, 'CORE MATERIAL', headerColor);
    
    y = sy - header_pad;
    xl = sx + 0.012; xv = sx + sw - 0.012;
    
    text(sx + sw/2, y, material.name, 'FontSize', 12, 'FontWeight', 'bold', ...
        'Color', accentColor, 'HorizontalAlignment', 'center'); y = y - dy;
    draw_line(xl, xv, y, 'k', sprintf('%.2e', material.k), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'α', sprintf('%.3f', material.alpha), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'β', sprintf('%.3f', material.beta), labelColor, valueColor);

    % =====================================================================
    % SECTION 4: WINDING (Top Right - Lower)
    % =====================================================================
    sy_winding = row1_top - h_material - col_gap;
    sh = h_winding;
    draw_section_box(sx, sy_winding - sh, sw, sh, boxColor);
    draw_header(sx + sw/2, sy_winding - 0.022, 'WINDING', headerColor);
    
    y = sy_winding - header_pad;
    draw_line(xl, xv, y, 'Stackup', stackup, labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'Copper', sprintf('%.0f / %.0f μm', t_cu_pri*1e6, t_cu_sec*1e6), labelColor, valueColor);

    % =====================================================================
    % SECTION 5: TRANSFORMER GEOMETRY (Bottom Left)
    % =====================================================================
    sx = col1_x; sy = row2_top; sw = col_w; sh = h_geometry;
    draw_section_box(sx, sy - sh, sw, sh, boxColor);
    draw_header(sx + sw/2, sy - 0.022, 'GEOMETRY', headerColor);
    
    y = sy - header_pad;
    xl = sx + 0.012; xv = sx + sw - 0.012;
    
    shape = opt.opt_design.centerpost_shape;
    text(sx + sw/2, y, sprintf('[ %s ]', shape), ...
        'FontSize', 9, 'Color', accentColor, 'HorizontalAlignment', 'center'); 
    y = y - dy;
    
    draw_line_sub(xl, xv, y, 'l', 'core', '', sprintf('%.2f mm', opt.opt_design.l_core*1e3), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'h', 'core', '', sprintf('%.2f mm', opt.opt_design.h_core*1e3), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'w', 'core', '', sprintf('%.2f mm', opt.opt_design.w_core*1e3), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'l', 'g', '', sprintf('%.3f mm', opt.opt_design.lg*1e3), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'Volume', sprintf('%.2f cm³', opt.opt_design.V_total*1e6), labelColor, valueColor); y = y - dy;
    draw_line(xl, xv, y, 'Footprint', sprintf('%.2f cm²', opt.opt_design.A_footprint*1e4), labelColor, valueColor);

    % =====================================================================
    % SECTION 6: TRANSFORMER LOSSES (Bottom Center)
    % =====================================================================
    sx = col2_x; sy = row2_top; sw = col_w; sh = h_losses;
    draw_section_box(sx, sy - sh, sw, sh, boxColor);
    draw_header(sx + sw/2, sy - 0.022, 'LOSSES', headerColor);
    
    y = sy - header_pad;
    xl = sx + 0.012; xv = sx + sw - 0.012;
    
    draw_line_sub(xl, xv, y, 'B', 'max', '', sprintf('%.4f T', opt.Bmax_opt), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'P', 'v', '', sprintf('%.1f kW/m³', opt.Pv_max_opt/1e3), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'P', 'core', '', sprintf('%.2f W', opt.P_core_min), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'P', 'cu', '', sprintf('%.2f W', opt.P_copper_min), labelColor, valueColor); y = y - dy;
    
    % Highlighted total
    y = y - 0.005;
    rect_h = dy + 0.01;
    rectangle('Position', [sx + 0.006, y - rect_h/2, sw - 0.012, rect_h], ...
        'FaceColor', [1, 0.95, 0.9], 'EdgeColor', accentColor, ...
        'LineWidth', 1.5, 'Curvature', 0.3);
    draw_line_sub(xl, xv, y, 'P', 'total', '', sprintf('%.2f W', opt.P_total_min), accentColor, accentColor);

    % =====================================================================
    % SECTION 7: THERMAL (Bottom Right - Upper)
    % =====================================================================
    sx = col3_x; sy = row2_top; sw = col_w; sh = h_thermal;
    draw_section_box(sx, sy - sh, sw, sh, boxColor);
    draw_header(sx + sw/2, sy - 0.022, 'THERMAL', headerColor);
    
    y = sy - header_pad;
    xl = sx + 0.012; xv = sx + sw - 0.012;
    
    draw_line_sub(xl, xv, y, 'R', 'dc,pri', '', sprintf('%.3f mΩ', opt.opt_design.Rdc_pri*1e3), labelColor, valueColor); y = y - dy;
    draw_line_sub(xl, xv, y, 'R', 'dc,sec', '', sprintf('%.3f mΩ', opt.opt_design.Rdc_sec*1e3), labelColor, valueColor); y = y - dy;
    
    if T_tx > 150
        tempColor = [0.8, 0.2, 0.2];
    elseif T_tx > 100
        tempColor = [0.9, 0.6, 0.1];
    else
        tempColor = [0.2, 0.6, 0.3];
    end
    draw_line_sub(xl, xv, y, 'T', 'xfmr', '', sprintf('%.1f °C', T_tx), labelColor, tempColor);

    % =====================================================================
    % EFFICIENCY BAR (Bottom Right - Lower)
    % =====================================================================
    eff = 100 * (1 - opt.P_total_min / Pmax);
    
    bar_margin = 0.008;
    bar_x = col3_x + bar_margin;
    bar_w = col_w - 2*bar_margin;
    bar_h = 0.045;
    bar_y = row2_top - h_thermal - col_gap - bar_h - 0.03;
    
    text(bar_x + bar_w/2, bar_y + bar_h + 0.022, 'EFFICIENCY', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', headerColor, ...
        'HorizontalAlignment', 'center');
    
    % Background
    rectangle('Position', [bar_x, bar_y, bar_w, bar_h], ...
        'FaceColor', [0.85, 0.85, 0.85], 'EdgeColor', [0.7, 0.7, 0.7], ...
        'LineWidth', 1, 'Curvature', 0.4);
    
    % Fill
    fill_w = bar_w * min(eff / 100, 1);
    if eff >= 99
        fillColor = [0.2, 0.7, 0.3];
    elseif eff >= 98
        fillColor = [0.5, 0.75, 0.25];
    else
        fillColor = [0.9, 0.6, 0.1];
    end
    rectangle('Position', [bar_x, bar_y, fill_w, bar_h], ...
        'FaceColor', fillColor, 'EdgeColor', 'none', 'Curvature', 0.4);
    
    text(bar_x + bar_w/2, bar_y + bar_h/2, sprintf('%.2f%%', eff), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

    % =====================================================================
    % FOOTER
    % =====================================================================
    text(0.5, 0.015, sprintf('Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM')), ...
        'FontSize', 8, 'Color', [0.5,0.5,0.5], 'HorizontalAlignment', 'center');

    hold off;
    axis([0 1 0 1]);
end

%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function draw_section_box(x, y, w, h, color)
    rectangle('Position', [x, y, w, h], ...
        'FaceColor', color, 'EdgeColor', [0.75, 0.75, 0.8], ...
        'LineWidth', 1, 'Curvature', 0.05);
end

function draw_header(x, y, txt, color)
    text(x, y, txt, 'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', color, 'HorizontalAlignment', 'center');
end

function draw_line(xl, xv, y, label, value, labelColor, valueColor)
    text(xl, y, label, 'FontSize', 9, 'Color', labelColor, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    text(xv, y, value, 'FontSize', 9, 'FontWeight', 'bold', 'Color', valueColor, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end

function draw_line_sub(xl, xv, y, base, sub, suffix, value, labelColor, valueColor)
    % Create label with subscript using TeX interpreter
    labelStr = sprintf('%s_{%s}%s', base, sub, suffix);
    text(xl, y, labelStr, 'FontSize', 9, 'Color', labelColor, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'tex');
    text(xv, y, value, 'FontSize', 9, 'FontWeight', 'bold', 'Color', valueColor, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end