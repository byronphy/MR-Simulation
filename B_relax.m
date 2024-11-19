function [M_end] = B_relax(M_start, t, M0, R1, R2)
% bloch_relax - compute the relaxation effect on the magnetization
%
% INPUTS
%	M_start: initial magnetization,3x1
%	t: duration [ms]
%	M0: equilibrium magnetization
%	R1: longitudinal relaxation time [kHz]
%	R2: transverse relaxation time [kHz]
%
% OUTPUTS
%   M_end - final magnetization

Arelax = [exp(-t*R2), 0, 0; ...
          0, exp(-t*R2), 0; ...
          0, 0, exp(-t*R1)];
brecover = [0; 0; M0*(1-exp(-t*R1))];
	
M_end = Arelax*M_start + brecover;

end