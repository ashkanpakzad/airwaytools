%% peripheral taper fitting
function [p_lumen, p_wall, fit_arc_length] = PeripheralTaper(bifurcation_fit, bifurcation_idx,logic_include, arc_length, lumen_area,wall_area)
% INPUT
%   bifurcation_fit = bifurcaiton number to fit peripheral taper rate from.
%   logic_include = vector of 0's and 1's where 0 indicates that the area
%   measurement is part of a bifurcation point.
%   arc_length = arc_length of each area measurement along airway.
%   lumen_area & wall_area = area measurements made along airway.
% OUTPUT
%   p_lumen % p_wall = fitted parameters of the peripheral airways. p(1)
%   represents the negative of the taper rate.
%   fit_arc_length = the arc_length values used for the fitting. provided
%   so that results can be visualised.

% EXAMPLE TO VISUALISE:
% line_lumen = p_lumen(1)*fit_arc_length+p_lumen(2);
% plot(fit_arc_length, line_lumen)

% fitting from the xth bifurcation (not including carina).
bifurcation_fit_idx = bifurcation_idx(bifurcation_fit);

% prepare copies of arrays from xth bifurcation
prefit_logic_include = logic_include(bifurcation_fit_idx:end);
prefit_arc_length = arc_length(bifurcation_fit_idx:end);
prefit_lumen_area = lumen_area(bifurcation_fit_idx:end);
prefit_wall_area = wall_area(bifurcation_fit_idx:end);

% identify points for peripheral fitting
fit_arc_length = prefit_arc_length(prefit_logic_include);
fit_lumen_area = prefit_lumen_area(prefit_logic_include);
fit_wall_area = prefit_wall_area(prefit_logic_include);

p_lumen = polyfit(fit_arc_length,fit_lumen_area,1);
p_wall = polyfit(fit_arc_length,fit_wall_area,1);
end