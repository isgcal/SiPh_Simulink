%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sen Lin <senlin@berkeley.edu>
% Integrated Systems Group, EECS, UC Berkeley
% 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
% open the library: SiPhotonics.slx
% open the testbench: ring_mod_link.slx
load_system('../SiPhotonics');
open_system('../ring_filter');

% link parameters
P_laser = 1;    % laser power
% model selection
opt = 1;      % default diode
% 1--linear, 2--sqrt, 3--diode

% technology
n0 = 2;         % effective index at f0
ng = 3;         % group index
a0 = 100;       % waveguide absorption coefficient (/m) 4.3dB/cm (intrinsic) 
Ne = 5e17;      % N-doping (/cm^3)
Nh = 5e17;      % P-doping (/cm^3)
Lj = 0.5e-4;    % junction geometry factor (cm)
a1 = -0.3;      % alpha 1st order coef. for linear/sqrt model
n1 = 8e-5;      % neff 1st order coef. for linear/sqrt model
Vbi = 0.8;      % built-in voltage, only for sqrt model
f0 = 230e12;    % measurement frequency for n0, a0

% device parameters microring modulator
L = 30e-6;          % ring perimeter
t_in = 0.99;        % input port coupling
t_drop = 0.995;      % drop port coupling

T0 = 300;   % measurement temp for parameters
Temp = 300; % current temp
nt = 2e-4;  % neff temperature coeff.
gamma_t = 1;% optical power to temp coeff. default = 0

% frequency sweep
fref = 229.96e12;  % simualtion reference frequency

fos = 0;

sim('ring_filter');
freq0 = fos0.Data;
len = length(freq0);
Tr = Tr0.Data(80:len);
lam = 3e8./(freq0(80:len)+fref)*1e9;
figure(1); hold on;
plot(lam, Tr);
xlabel('wavelength (nm)');