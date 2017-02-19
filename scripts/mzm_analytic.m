%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sen Lin <senlin@berkeley.edu>
% Integrated Systems Group, EECS, UC Berkeley
% 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ring_mod_analytic.m
clear all; close all;
% model selection
opt = 3;      % default diode
% 1--linear, 2--sqrt, 3--diode

% technology
n0 = 2;         % effective index at f0
ng = 3;         % group index
a0 = 100;       % absorption coefficient (/m) 4.3dB/cm (intrinsic)
Ne = 5e17;      % N-doping (/cm^3)
Nh = 5e17;      % P-doping (/cm^3)
Lj = 0.5e-4;    % junction geometry factor (cm)
a1 = -0.3;      % alpha 1st order coef. for linear/sqrt model
n1 = 1e-4;      % neff 1st order coef. for linear/sqrt model
Vbi = 0.8;      % built-in voltage, only for sqrt model
f0 = 230e12;    % measurement frequency for n0, a0
wf = 0.9;        % waveguide factor

% device parameters for MZM
L = 1e-3;             % arm length in meter
phase_os0 = 0.3;
phase_os1 = phase_os0-0.5*pi;   % quarature point

% frequency sweep
fref = 230.5e12;  % simualtion reference frequency
fstart = fref;
fend = 231e12;
fstep = 100e9;
f_swp = fstart:fstep:fend;
lambda = 3e8./f_swp;

mzm0 = MZModulator(opt, L, n0, ng, a0, Ne, Nh, Lj, wf, a1, n1, Vbi,f0);

Vdrive_0 = 0.5;
Vdrive_1 = -1.5;


% find the optimal phase_os0 to maximize OMA
laser_opt = 1300e-9;
flaser = 3e8/laser_opt; % used in simulink
phase_os0_swp = 0:pi/100:pi/2;
OMA_swp = zeros(1, length(phase_os0_swp));
j = 0;
for phase_os0 = phase_os0_swp
    j = j + 1;
    phase_os1 = phase_os0-0.6*pi;   % quarature point
    [~, OMA_swp(j), P1, P0] = mzm0.modulation(laser_opt, phase_os0, phase_os1, Vdrive_0, Vdrive_1);
end
[OMA_max,ind_max] = max(OMA_swp);
phase_os0 = phase_os0_swp(ind_max);
phase_os1 = phase_os0-0.5*pi;   % quarature point
[~, ~, P1, P0] = mzm0.modulation(laser_opt, phase_os0, phase_os1, Vdrive_0, Vdrive_1);

fprintf('-------MZ Modulator Analytical Model-------\n');
fprintf('Optimal laser wavelength:  %4.4f nm\n', 1e9*laser_opt);
fprintf('Optimal laser frequency:  %3.5f THz\n', 1e-12*3e8/laser_opt);
fprintf('Max. OMA: %1.3f \n', OMA_max);
fprintf('Phase_os0 = %1.3f, Phase_os1 = %1.3f \n', phase_os0, phase_os1);
fprintf('Bit-1 power P1: %1.3f \n', P1);
fprintf('Bit-0 power P0: %1.3f \n', P0);

[~,~,Pt0, Pd0] = mzm0.tf(lambda, phase_os0, phase_os1, Vdrive_0, Vdrive_1);
[~,~,Pt1, Pd1] = mzm0.tf(lambda, phase_os0, phase_os1, Vdrive_1, Vdrive_0);

close all;
figure(1); 
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 800, 400]);hold on; 
hold on; grid on;
set(gca, 'FontSize', 18, 'LineWidth', 2); 
set(gca,'Box','on');
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
plot(lambda*1e9, Pt0, 'r','LineWidth',3); 
plot(lambda*1e9, Pt1, 'b','LineWidth',3); 
axis([min(lambda*1e9) max(lambda*1e9) 0 1]);
xlabel('Wavelength (nm)');
ylabel('Power Transmission');
legend('bit0', 'bit1','Location', 'northwest');

plot([1e9*laser_opt, 1e9*laser_opt],[0, 1],'m-.','LineWidth',2);

% find the neff and alpha under different biases
[neff0, alpha0] = mzm0.material(laser_opt, Vdrive_0);
[neff1, alpha1] = mzm0.material(laser_opt, Vdrive_1);
fprintf('Vdrive_0 = %g V, ', Vdrive_0);
fprintf('alpha0 = %g/m, ', alpha0);
fprintf('neff0 = %g\n', neff0);
fprintf('Vdrive_1 = %g V, ', Vdrive_1);
fprintf('alpha1 = %g/m, ', alpha1);
fprintf('neff1 = %g\n', neff1);
fprintf('loss = %g dB/cm\n', 10*log10(exp(-(alpha1+alpha0)/2e2)));
fprintf('VpiL = %g Vcm\n', 1e2*0.5*laser_opt/((neff1-neff0)/(Vdrive_0-Vdrive_1)));







