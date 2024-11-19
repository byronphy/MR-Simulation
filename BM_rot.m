function [M_end] = BM_rot(M_start, dt, B)
% Compute the rotation of the water proton magnetization
%
% INPUTS
%	M_start: initial magnetization, 3x1 for water proton pool
%	dt: duration [ms]
%	B = [Bx, By, Bz] - magnetic field [mT]
%
% OUTPUTS
%   M_end - final magnetization

GAMMA = 42.58; % [kHz/mT] for proton
precess_angle = -2*pi*GAMMA * norm(B) * dt; % [rad], minus for the MR precessing
theta = acos(B(3) / (norm(B)+eps));
phi = atan2(B(2), B(1));

M_end = Rz(phi)*Ry(theta)*Rz(precess_angle)*Ry(-theta)*Rz(-phi)* M_start;

end