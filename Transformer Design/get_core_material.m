% Get Core Material Parameters
%
% Author:  Julia Estrin
% Date:    02-27-2026
% Updated: 03-02-2026 - full steinmetz params 
%
% Description:
%   Returns Steinmetz parameters and relative permeability for common
%   ferrite core materials. Supports easy material selection in design scripts.
%
% Inputs:
%   material_name - String identifying the material
%                   Options: 'ML91S', 'F80', '3F4'
%
% Outputs:
%   params - Struct containing:
%       .k    - Steinmetz coefficient k  [W·s^beta/(T^beta·m^3)]
%       .beta - Steinmetz exponent beta  [-]
%       .alpha - Steinmetz exponent beta  [-]
%       .uc   - Relative permeability    [-]
%       .name - Material name string
%
% Usage Example:
%   material = get_core_material('ML91S');
%   k    = material.k;
%   beta = material.beta;
%   alpha = material.alpha;
%   uc   = material.uc;

function params = get_core_material(material_name)

    switch upper(material_name)
        
        case {'ML91S', 'PROTERIALS_ML91S', 'PROTERIALS ML91S'}
            % Proterials ML91S - High frequency ferrite
            % Optimal for 500 kHz - 3 MHz
            params.k    = 1.12893E-08;
            params.beta = 3.095989317;
            params.alpha = 2.8298278;
            params.uc   = 900;
            params.name = 'Proterials ML91S';
        
        case {'F80', 'FERRITE_80', 'FERRITE 80'}
            % Ferrite 80 - General purpose high frequency ferrite
            % Optimal for 200 kHz - 2 MHz
            params.k    = 6.43007E-05;
            params.beta = 2.549625682;
            params.alpha = 2.15651264669503;
            params.uc   = 600;
            params.Bsat = 0.3; % [T]
            params.name = 'Fairite 80';
        
              
        % case {'3F46', 'FERROXCUBE_3F46', 'FERROXCUBE 3F46'}
        %     params.k    = 8.69811E+07;
        %     params.beta = 2.191149652;
        %     params.uc   = 750;
        %     params.name = 'Ferroxcube 3F46';
        
        otherwise
            error('get_core_material:UnknownMaterial', ...
                ['Unknown material: %s\n' ...
                 'Available materials: ML91S, F80, 3F4'], ...
                material_name);
    end

end