%% Rotating Frame, Bloch Equation
clc; clear; close all;

GAMMA = 42.58; % [MHz/T, kHz/mT]
B0 = 7; % [T]
M0 = 1;
M_start = [0,0,M0].';
R1 = 1/2000; % [kHz]
R2 = 1/80; %[kHz]

dt = 0.01; % [ms]
% 
% % Hard Pulse
% flipAngle = 10; % [deg]
% pulse_dur = 0.100; % [ms]
% pulse_delay = 0.200; % [ms]
% Nt_dur = ceil(pulse_dur/dt);
% Nt_delay = ceil(pulse_delay/dt);
% N_unit = 30;
% pulse_unit = [ones(1,Nt_dur),zeros(1,Nt_delay),-ones(1,Nt_dur),zeros(1,Nt_delay)];
% pulse_all = repmat(pulse_unit, 1, N_unit);
% pulse_all = [-0.5*ones(1,Nt_dur),zeros(1,Nt_delay),pulse_all,0.5*ones(1,Nt_dur),zeros(1,Nt_delay)];
% B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*pulse_dur)*1e3; % [mT]
% B1p_all = pulse_all*B1_scale; % [mT]
% t_end = 2*(pulse_dur+pulse_delay)*(N_unit+1); % [ms]

% % test Hard Pulse
% flipAngle = 90; % [deg]
% pulse_dur = 3; % [ms]
% Nt_dur = ceil(pulse_dur/dt);
% B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*pulse_dur); % [mT]
% pulse_all = [ones(1,Nt_dur),zeros(1,Nt_dur)];
% B1p_all = pulse_all*B1_scale; % [mT]
% t_end = dt*length(B1p_all); % [ms]

% test Gauss pulse
load('GAUSS5120_B375.mat');
B1p_gauss = rf.waveform;
B1p_gauss = B1p_gauss(:,1).*(cos(B1p_gauss(:,2))+1j*sin(B1p_gauss(:,2)));
flipAngle = 30; % [deg]
B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*sum(B1p_gauss*dt)); % [mT]
B1p_gauss = B1p_gauss.'*B1_scale;

pulse_delay = 100; % [ms]
Nt_delay = ceil(pulse_delay/dt);
N_unit = 2;
B1p_unit = [B1p_gauss, zeros(1,Nt_delay)];
B1p_all = repmat(B1p_unit, 1, N_unit);
t_end = dt*length(B1p_all);



% time domain
t_begin = 0;
t_all = t_begin:dt:t_end-dt; % [ms]
Nt = length(t_all); % 

BW_ppm = 5; % [ppm]
BW = BW_ppm*GAMMA*B0*1e-3; % [kHz]
%BW = 2; % [kHz] off-resonance frequency bandwidth (e.g. the linear gradient field)
df = linspace(-BW,BW); % [kHz]

M = repmat(M_start, [1, length(df)]);
for I_t = 1:Nt
    for I_f = 1:length(df)
        Bz = df(I_f)/GAMMA; % [mT] linear Bz field
        M(:,I_f) = BM_rot(M(:,I_f), dt/2, [real(B1p_all(I_t)),imag(B1p_all(I_t)),Bz]);
        M(:,I_f) = B_relax(M(:,I_f), dt, M0, R1, R2);
        M(:,I_f) = BM_rot(M(:,I_f), dt/2, [real(B1p_all(I_t)),imag(B1p_all(I_t)),Bz]);
    end
end

figure;
subplot(131)
plot(t_all,B1p_all)
xlabel('time (ms)'), ylabel('B_1^+ (mT)')
title('RF Pulse in Rotating Frame')
subplot(132)
plot(df, M(3,:))
title('Frequency profile')
xlabel('Off-resonance frequency (kHz)'), legend('M_Z')
subplot(133)
plot(df,sqrt(M(1,:).^2 + M(2,:).^2))
title('Frequency profile')
xlabel('Off-resonance frequency (kHz)'), legend('|M_{XY}|')



%%
%% MT time evolution
%clc; clear; close all;

% GAMMA = 42.58; % [MHz/T, kHz/mT]
% M0 = 1;
% M_start = [0;0;M0]; % initial magnetization vecotr [Mx,My,Mz]
% 
% R1 = 1/2000; % [kHz]
% R2 = 1/80; %[kHz]
% 
% % Hard Pulse
% dt = 0.01; % [ms]
% 
% % % Rectangular pulses (rotating frame)
% % flipAngle = 10; % [deg]
% % pulse_dur = 1; % [ms]
% % pulse_delay = 1; % [ms]
% % Nt_dur = ceil(pulse_dur/dt);
% % Nt_delay = ceil(pulse_delay/dt);
% % N_unit = 2;
% % pulse_unit = [ones(1,Nt_dur),zeros(1,Nt_delay),-ones(1,Nt_dur),zeros(1,Nt_delay)];
% % pulse_all = repmat(pulse_unit, 1, N_unit);
% % pulse_all = [-0.5*ones(1,Nt_dur),zeros(1,Nt_delay),pulse_all,0.5*ones(1,Nt_dur),zeros(1,Nt_delay)];
% % B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*pulse_dur)*1e3; % [mT]
% % B1p_all = pulse_all*B1_scale; % [mT]
% % t_end = 2*(pulse_dur+pulse_delay)*(N_unit+1); % [ms]
% 
% % test Hard Pulse
% flipAngle = 30; % [deg]
% pulse_dur = 3; % [ms]
% Nt_dur = ceil(pulse_dur/dt);
% B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*pulse_dur); % [mT]
% pulse_all = [ones(1,Nt_dur),zeros(1,Nt_dur)];
% B1p_all = pulse_all*B1_scale; % [mT]
% t_end = dt*length(B1p_all); % [ms]

% % time domain
% t_begin = 0;
% t_all = t_begin:dt:t_end-dt; % [ms]
% Nt = length(t_all); % 

M = repmat(M_start, [1,Nt]);
for I_t = 1:Nt-1
    M(:,I_t+1) = BM_rot(M(:,I_t), dt/2, [real(B1p_all(I_t)),imag(B1p_all(I_t)),0]);
    M(:,I_t+1) = B_relax(M(:,I_t+1), dt, M0, R1, R2);
    M(:,I_t+1) = BM_rot(M(:,I_t+1), dt/2, [real(B1p_all(I_t)),imag(B1p_all(I_t)),0]);
end

figure;
subplot(311)
plot(t_all,B1p_all)
xlabel('time (ms)'), ylabel('B_1^+ (mT)')
title('RF Pulse in Rotating Frame')
subplot(312)
plot(t_all, M(3,:)),hold on;
title('Time profile')
xlabel('time (ms)'), legend('M_Z')
subplot(313)
plot(t_all,sqrt(M(1,:).^2 + M(2,:).^2))
title('Time profile')
xlabel('time (ms)'), legend('|M_{XY}|')
