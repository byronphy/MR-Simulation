
B0 = 7; % [T]
GAMMA = 42.58; % [MHz/T], [kHz/mT]

BW_ppm = 5; % [ppm]
BW = BW_ppm*GAMMA*B0*1e-3; % [kHz]
%BW = 2; % [kHz] off-resonance frequency bandwidth (e.g. the linear gradient field)
df = linspace(-BW,BW); % [kHz]
Nf = 1000;
df = (-BW:2*BW/Nf:BW);

flipAngle = 30; % [deg]
T_rf = 5.12; % [ms]
B1 = (flipAngle*pi/180)/(2*pi*GAMMA*T_rf); % [mT]

T_d = 5; % [ms]

theta = acos(df./sqrt(GAMMA^2*B1^2+df.^2));
c = cos(theta);
s = sin(theta);

f_df = T_d*df*2*pi;
cf = cos(f_df);
sf = sin(f_df);

g_df = T_rf*sqrt(GAMMA^2*B1^2+df.^2)*2*pi;
cg = cos(g_df);
sg = sin(g_df);

gh_df = g_df/2;
cgh = cos(gh_df);
sgh = sin(gh_df);

Mz = (c.^2+cg.*s.^2).^2-s.^2.*sg.*(-2*c.*sf.*sgh.^2+cf.*sg)-c.*(-1+cg).*s.^2.*(2*cf.*c.*sgh.^2+sf.*sg);
figure;
plot(df,Mz);
title('Theoretical Frequency profile')
xlabel('Off-resonance frequency (kHz)'), legend('M_Z')