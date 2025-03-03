function [r_fin, r_d] = compartment_radius(stretch, L_t, r_ini)
% FUNCTION: Computes the final radius of a compartment based on stretch.
% INPUTS:
%   - stretch : Stretch factor (unitless, should be >1 for elongation)
%   - L_t     : Tissue length undeformed (unit: cm)
%   - r_ini   : Initial compartment radius (unit: cm)
% OUTPUTS:
%   - r_fin   : Final compartment radius after stretch (unit: cm)
%   - r_d     : Reduction in radius due to stretch (unit: cm)


% -------------------- Stretch Equilibrium Adjustment --------------------

stretch_eq_num = 0.005;  % Adjustment factor for stretch equilibrium

% -------------------- Compute Radius Reduction Due to Stretch --------------------

if stretch > 1  % Stretch must be greater than 1 for elongation

    % Compute the deformed length (elongation of the tissue)
    L_d = (stretch - stretch_eq_num) * L_t;

    % Ensure L_d does not shrink below original tissue length
    if L_d <= L_t
        L_d = L_t;
    end

    % Compute omega values for geometry-based radius estimation
    omega_1 = L_d / 4;
    omega_d = (omega_1 - (L_t / 4)) / 2;
    omega_2 = sqrt(power(omega_1, 2) - power(omega_d, 2));
    omega_3 = 2 * omega_2;

    % Compute radius reduction due to elongation
    r_d = sqrt(power(omega_3, 2) - power(L_t / 2, 2));

    % Compute the final radius after deformation
    r_fin = r_ini - r_d;

else
    % No stretching occurs; maintain original radius
    r_fin = r_ini;
    r_d = 0;
end

end