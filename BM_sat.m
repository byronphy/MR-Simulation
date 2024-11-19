function [M_end] = BM_sat(M_start, dt, R_RFB)
% Compute the saturation of the macromolecular proton magnetization
%
% INPUTS
%	M_start: initial magnetization, 1x1 for macromolecular proton pool
%	dt: duration [ms]
%	R_RFB: saturation rate [kHz]
%
% OUTPUTS
%   M_end - final magnetization

M_end = exp(-R_RFB*dt)* M_start;

end