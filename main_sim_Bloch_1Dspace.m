%% Rotating Frame, Bloch Equation
clc; clear; close all;

GAMMA = 42.58; % [MHz/T, kHz/mT]
B0 = 7; % [T]
M0 = 1;
M_start = [0,0,M0].';
R1 = 1/2000; % [kHz]
R2 = 1/80; %[kHz]

dt = 0.01; % [ms]

% test Hard Pulse
flipAngle = 30; % [deg]
pulse_dur = 5.12; % [ms]
Nt_dur = ceil(pulse_dur/dt);
B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*pulse_dur); % [mT]
B1p_rect = ones(1,Nt_dur)*B1_scale; % [mT]

% % test Gauss pulse
% load('GAUSS5120_B375.mat');
% B1p_gauss = rf.waveform;
% B1p_gauss = B1p_gauss(:,1).*(cos(B1p_gauss(:,2))+1j*sin(B1p_gauss(:,2)));
% flipAngle = 30; % [deg]
% B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*sum(B1p_gauss*dt)); % [mT]
% B1p_gauss = B1p_gauss.'*B1_scale;

pulse_delay = 5; % [ms]
Nt_delay = ceil(pulse_delay/dt);


% 1D Space
radius = 2; % [mm]
N_x = 9; % number of magnetizations, int
dx = linspace(-radius,radius,N_x); % [mm]
dx = dx(1 : (end - 1)); % Asymmetry
% Crusher Gradient
crush_dur = 3; % [ms]
Nt_cru = ceil(crush_dur/dt);
crush_str = 10*1e-3; % [mT/mm]
dBz_x = crush_str*dx; % [mT]


% unit
N_unit = 2;
%B1p_unit = [B1p_gauss, zeros(1,Nt_delay)];
B1p_unit = [B1p_rect, zeros(1,Nt_delay), zeros(1,Nt_cru)];
B1p_all = repmat(B1p_unit, 1, N_unit);

dBz_x_unit = repmat(zeros(1,length(dx)),[Nt_dur+Nt_delay+Nt_cru,1]); % [mT]
dBz_x_unit(ceil(Nt_dur+Nt_delay+1):end,:)=repmat(dBz_x,[Nt_cru,1]); % [mT]
dBz_x_all = repmat(dBz_x_unit, N_unit,1);

dGx_unit = zeros(1,length(dBz_x_unit));
dGx_unit(ceil(Nt_dur+Nt_delay+1):end)=crush_str; % [mT/mm]
dGx_all = repmat(dGx_unit, 1,N_unit);


% time domain
t_begin = 0;
t_end = dt*length(B1p_all);
t_all = t_begin:dt:t_end-dt; % [ms]
Nt = length(t_all); % 

% frequency domain
BW_ppm = 5; % [ppm]
BW = BW_ppm*GAMMA*B0*1e-3; % [kHz]
%BW = 2; % [kHz] off-resonance frequency bandwidth (e.g. the linear gradient field)
df = linspace(-BW,BW); % [kHz]
Nf = 500;
df = (-BW:2*BW/Nf:BW);

M = repmat(M_start, [1, length(df), length(dx), Nt]);
for I_t = 1:Nt-1
    for I_f = 1:length(df)
        for I_x = 1:length(dx)
            Bz = df(I_f)/GAMMA; % [mT] linear Bz field
            %M(:,I_f) = BM_rot(M(:,I_f), dt/2, [real(B1p_all(I_t)),imag(B1p_all(I_t)),Bz]);
            %M(:,I_f) = B_relax(M(:,I_f), dt, M0, R1, R2);
            %M(:,I_f) = BM_rot(M(:,I_f), dt/2, [real(B1p_all(I_t)),imag(B1p_all(I_t)),Bz]);
            M(:,I_f,I_x,I_t+1) = BM_rot(M(:,I_f,I_x,I_t), dt, [real(B1p_all(I_t)),imag(B1p_all(I_t)),Bz+dBz_x_all(I_t,I_x)]);
        end
    end
end

% Sum along the 1D space
M_voxel(:,:) = sum(M(:,:,:,end),3)/length(dx);

M_resonance = squeeze(M(:,Nf/2+1,:,:));
M_time(:,:) = sum(M_resonance,2)/length(dx);

%
figure(1);
subplot(131)
plot([t_all(1)-eps, t_all, t_all(end)+eps],[0,B1p_all,0]);hold on;
plot([t_all(1)-eps, t_all, t_all(end)+eps],[0,dGx_all/10,0]);hold off;
xlabel('time (ms)'), legend({'B_1^+ (mT)','Gx (10 mT/mm)'})
title('RF Pulse in Rotating Frame')
subplot(132)
plot(df, M_voxel(3,:))
title('Frequency profile')
xlabel('Off-resonance frequency (kHz)'), legend('M_Z')
subplot(133)
plot(df,sqrt(M_voxel(1,:).^2 + M_voxel(2,:).^2))
title('Frequency profile')
xlabel('Off-resonance frequency (kHz)'), legend('|M_{XY}|')

figure(2);
subplot(311)
plot([t_all(1)-eps, t_all, t_all(end)+eps],[0,B1p_all,0]);hold on;
plot([t_all(1)-eps, t_all, t_all(end)+eps],[0,dGx_all/10,0]);hold off;
xlabel('time (ms)'), legend({'B_1^+ (mT)','Gx (10 mT/mm)'})
title('RF Pulse in Rotating Frame')
subplot(312)
plot([t_all(1)-eps, t_all, t_all(end)+eps], [0,M_time(3,:),0]);
title('Time profile')
xlabel('time (ms)'), legend('M_Z')
subplot(313)
plot([t_all(1)-eps, t_all, t_all(end)+eps],[0,sqrt(M_time(1,:).^2 + M_time(2,:).^2),0])
title('Time profile')
xlabel('time (ms)'), legend('|M_{XY}|')
