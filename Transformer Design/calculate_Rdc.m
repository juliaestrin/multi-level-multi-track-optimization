% Calculate Equivalent DC Resistance of a Planar Winding (Round or Rectangular)
%
% Author:  Julia Estrin
% Date:    03-27-2026
%
% Description:
%   Wrapper function that computes the equivalent DC resistance of a planar
%   winding, automatically selecting the appropriate geometry model based on
%   the winding_shape parameter. Supports both rectangular and round (annular)
%   winding geometries.
%
% Inputs:
%   n             - Index of the outermost turn (loop runs from turn 0 to turn n)
%   sig_cu        - Electrical conductivity of copper [S/m]
%   t_cu          - Copper trace thickness [m]
%   winding_shape - Winding geometry: 'round' or 'square'/'rectangular'
%   varargin      - Shape-dependent geometry parameters:
%
%                   For 'square' or 'rectangular':
%                     a_0   - Inner winding length in x-direction [m]
%                     a_n   - Outer winding length in x-direction [m]
%                     b_0   - Inner winding length in y-direction [m]
%                     b_n   - Outer winding length in y-direction [m]
%
%                   For 'round':
%                     r_0   - Inner winding radius [m]
%                     r_n   - Outer winding radius [m]
%
% Outputs:
%   Rdc           - Equivalent DC resistance of the winding [Ohm]
%
% Usage Examples:
%   % Rectangular winding
%   Rdc = calculate_Rdc(100, 5.8e7, 35e-6, 'square', 0.01, 0.05, 0.01, 0.05);
%
%   % Round winding
%   Rdc = calculate_Rdc(100, 5.8e7, 35e-6, 'round', 5e-3, 25e-3);

function Rdc = calculate_Rdc(n, sig_cu, t_cu, winding_shape, varargin)

    % --- Input validation ---
    winding_shape = lower(winding_shape);
    
    switch winding_shape
        case {'square', 'rectangular', 'rect'}
            % Rectangular winding requires 4 geometry parameters
            if length(varargin) < 4
                error('Rectangular winding requires 4 geometry parameters: a_0, a_n, b_0, b_n');
            end
            a_0 = varargin{1};
            a_n = varargin{2};
            b_0 = varargin{3};
            b_n = varargin{4};
            
            Rdc = calculate_Rdc_rectangle(n, a_0, a_n, b_0, b_n, sig_cu, t_cu);
            
        case {'round', 'circular', 'annular'}
            % Round winding requires 2 geometry parameters
            if length(varargin) < 2
                error('Round winding requires 2 geometry parameters: r_0, r_n');
            end
            r_0 = varargin{1};
            r_n = varargin{2};
            
            Rdc = calculate_Rdc_round(n, r_0, r_n, sig_cu, t_cu);
            
        otherwise
            error('winding_shape must be ''round'' or ''square''');
    end

end