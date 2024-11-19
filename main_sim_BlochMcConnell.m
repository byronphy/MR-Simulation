%% MT spectrum
clc; clear; close all;

GAMMA = 42.58; % [MHz/T, kHz/mT]
fraction = 0.1;
M0w = 1-fraction;
M0m = fraction;
M_start = [0;0;M0w;M0m]; % initial magnetization vecotr [Mxw,Myw,Mzw,Mzm]

K = 0.00145; % [kHz], 0.00145
K = 0;
R1w = 1/2500; % [kHz]
R2w = 1/50; %[kHz]
Kw = K/(1-fraction); %[kHz]

R1m = 1/500; % [kHz]
R2m = 1/0.07; % [kHz]
Km = K/fraction; %[kHz]

% Hard Pulse
dt = 0.0010; % [ms]

% Rectangular pulses (rotating frame)
flipAngle = 10;  % [deg]
pulse_dur = 0.100; % [ms]
pulse_delay = 0.200; % [ms]
Nt_dur = ceil(pulse_dur/dt);
Nt_delay = ceil(pulse_delay/dt);
N_unit = 30;
pulse_unit = [ones(1,Nt_dur),zeros(1,Nt_delay),-ones(1,Nt_dur),zeros(1,Nt_delay)];
pulse_all = repmat(pulse_unit, 1, N_unit);
pulse_all = [-0.5*ones(1,Nt_dur),zeros(1,Nt_delay),pulse_all,0.5*ones(1,Nt_dur),zeros(1,Nt_delay)];
B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*pulse_dur); % [mT]
B1p_all = pulse_all*B1_scale; % [mT]
t_end = 2*(pulse_dur+pulse_delay)*(N_unit+1); % [ms]

% % test Hard Pulse
% flipAngle = 90; % [deg]
% pulse_dur = 0.1; % [ms]
% Nt_dur = ceil(pulse_dur/dt);
% B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*pulse_dur); % [mT]
% pulse_all = [ones(1,Nt_dur),zeros(1,Nt_dur)];
% B1p_all = pulse_all*B1_scale; % [mT]
% t_end = dt*length(B1p_all); % [ms]

% time domain
t_begin = 0;
t_all = t_begin:dt:t_end-dt; % [ms]
Nt = length(t_all); % 

BW = 5; % [kHz] off-resonance frequency bandwidth (e.g. the linear gradient field)
df = linspace(-BW,BW); % [kHz]

M = repmat(M_start, [1, length(df)]);
for I_t = 1:Nt
    for I_f = 1:length(df)
        Bz = df(I_f)/GAMMA; % [mT]
        R_RFB = (B1p_all(I_t)*2*pi*GAMMA)^2/R2m; % [kHz, rad]
        M(1:3,I_f) = BM_rot(M(1:3,I_f), dt/2, [B1p_all(I_t),0,Bz]);
        M(4,I_f) = BM_sat(M(4,I_f), dt/2, R_RFB);
        M(:,I_f) = BM_relaxMT(M(:,I_f), dt, M0w, M0m, R2w, R1w, R1m, Kw, Km);
        M(4,I_f) = BM_sat(M(4,I_f), dt/2, R_RFB);
        M(1:3,I_f) = BM_rot(M(1:3,I_f), dt/2, [B1p_all(I_t),0,Bz]);
    end
end

figure;
subplot(131)
plot([t_all(1)-eps,t_all,t_all(end)+eps],[0,B1p_all,0])
xlabel('time (ms)'), ylabel('B_1^+ (mT)')
title('Hard RF Pulse in Rotating Frame')
subplot(132)
plot(df, M(3,:)),hold on;
plot(df, M(4,:)),hold off;
title('Frequency profile')
xlabel('Off-resonance frequency (kHz)'), ylabel('M_z')
legend({'M_z^w','M_z^m'}),legend boxoff;
subplot(133)
plot(df,M(1,:)),hold on;
plot(df,M(2,:)),hold on;
title('Frequency profile')
xlabel('Off-resonance frequency (kHz)'), ylabel('M_{x,y}')
legend({'M_x^w','M_y^w'}),legend boxoff;


%% MT time evolution
%clc; clear; close all;

% GAMMA = 42.58; % [MHz/T, kHz/mT]
% fraction = 0.1;
% M0w = 1-fraction;
% M0m = fraction;
% M_start = [0;0;M0w;M0m]; % initial magnetization vecotr [Mxw,Myw,Mzw,Mzm]
% 
% K = 0.00145; % [kHz], 0.00145
% %K = 0;
% R1w = 1/2000; % [kHz]
% R2w = 1/80; %[kHz]
% Kw = K/(1-fraction); %[kHz]
% 
% R1m = 1/500; % [kHz]
% R2m = 1/0.07; % [kHz]
% Km = K/fraction; %[kHz]
% 
% % Hard Pulse
% dt = 0.01; % [ms]

% % Rectangular pulses (rotating frame)
% flipAngle = 10; % [deg]
% pulse_dur = 0.100; % [ms]
% pulse_delay = 0.200; % [ms]
% Nt_dur = ceil(pulse_dur/dt);
% Nt_delay = ceil(pulse_delay/dt);
% N_unit = 30;
% pulse_unit = [ones(1,Nt_dur),zeros(1,Nt_delay),-ones(1,Nt_dur),zeros(1,Nt_delay)];
% pulse_all = repmat(pulse_unit, 1, N_unit);
% pulse_all = [-0.5*ones(1,Nt_dur),zeros(1,Nt_delay),pulse_all,0.5*ones(1,Nt_dur),zeros(1,Nt_delay)];
% B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*pulse_dur); % [mT]
% B1p_all = pulse_all*B1_scale; % [mT]
% t_end = 2*(pulse_dur+pulse_delay)*(N_unit+1); % [ms]

% % test Hard Pulse
% flipAngle = 90; % [deg]
% pulse_dur = 100; % [ms]
% Nt_dur = ceil(pulse_dur/dt);
% B1_scale = (flipAngle*pi/180)/(2*pi*GAMMA*pulse_dur); % [mT]
% pulse_all = [ones(1,Nt_dur),zeros(1,Nt_dur)];
% B1p_all = pulse_all*B1_scale; % [mT]
% t_end = dt*length(B1p_all); % [ms]

% time domain
t_begin = 0;
t_all = t_begin:dt:t_end-dt; % [ms]
Nt = length(t_all); % 

M = repmat(M_start, [1,Nt]);
for I_t = 1:Nt-1
    Bz = 0; % [mT]
    %Bz = 1.67/GAMMA; % [mT]
    R_RFB = (B1p_all(I_t)*GAMMA)^2/R2m; % [kHz, rad]
    M(1:3,I_t+1) = BM_rot(M(1:3,I_t), dt/2, [B1p_all(I_t),0,Bz]);
    M(4,I_t+1) = BM_sat(M(4,I_t+1), dt/2, R_RFB);
    M(:,I_t+1) = BM_relaxMT(M(:,I_t+1), dt, M0w, M0m, R2w, R1w, R1m, Kw, Km);
    M(4,I_t+1) = BM_sat(M(4,I_t+1), dt/2, R_RFB);
    M(1:3,I_t+1) = BM_rot(M(1:3,I_t+1), dt/2, [B1p_all(I_t),0,Bz]);
end

figure;
subplot(311)
plot(t_all,B1p_all)
xlabel('time (ms)'), ylabel('B_1^+ (mT)')
title('Hard RF Pulse in Rotating Frame')
subplot(312)
plot(t_all, M(3,:)),hold on;
plot(t_all, M(4,:)),hold off;
title('Time profile')
xlabel('time (ms)'), legend('M_Z'),legend boxoff;
subplot(313)
plot(t_all,M(1,:)),hold on;
plot(t_all,M(2,:)),hold on;
title('Time profile')
xlabel('time (ms)'), legend({'M_x','M_y'}),legend boxoff;
