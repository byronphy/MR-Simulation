function [M_end] = BM_relaxMT(M_start, dt, M0w, M0m, R2w, R1w, R1m, Kw, Km)
% Compute the relaxation and MT seffect on the magnetization
%
% INPUTS
%	M_start: initial magnetization, 4x1
%	dt: duration [ms]
%	M0w: equilibrium magnetization for water proton
%	M0w: equilibrium magnetization for macromolecular proton
%	R2w: transverse relaxation rate for water proton [kHz]
%	R1w: longitudinal relaxation rate for water proton [kHz]
%	R1m: longitudinal relaxation rate for macromolecular proton [kHz]
%	Kw: exchange rate for water -> macromolecular [kHz]
%	Km: exchange rate for macromolecular -> water [kHz]
%
% OUTPUTS
%   M_end - final magnetization, 4x1

% Splitting method for relaxation and MT effect
M_end(1) = exp(-R2w*dt)*M_start(1);
M_end(2) = exp(-R2w*dt)*M_start(2);

A_mt = [-R1w-Kw, Km;...
        Kw, -R1m-Km];
b_mt = [M0w*R1w; M0m*R1m];
if det(A_mt) == 0
    mb = [0,0,0,0];
else
    mb = A_mt \ b_mt;  % A_mt * mb=b_mt
end

[V, D] = eig(A_mt);
V_det = V(1,1)*V(2,2)-V(1,2)*V(2,1);
V_inv = [V(2,2),-V(1,2);...
         -V(2,1),V(1,1)]/V_det;
exp_D_dt = [exp(D(1,1)*dt),0;...
            0,exp(D(2,2)*dt)];
M_end(3:4) = V*exp_D_dt*V_inv*(M_start(3:4)+mb)-mb;

end