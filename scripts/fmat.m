%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sen Lin <senlin@berkeley.edu>
% Integrated Systems Group, EECS, UC Berkeley
% 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [neff, alpha] = fmat(fref, V, val, c0)
% material properties

model_select = val(1); 
% model selection
% 1--linear, 2--sqrt, 3--diode
n0 = val(2);  % effective refractive index measured at f0
ng = val(3);  % group index
a0 = val(4);  % absorption coefficient in m^-1
n1 = val(5);  % positive
a1 = val(6);  % negative
Nd = val(7);  % n-doping
Na = val(8);  % p-doping
Vbi = val(9); % build-in voltage
Lj = val(10); % geometry/mode factor (cm)
f0 = val(11); % measurement frequency
wf = val(12);  % waveguide factor

% define constants
eps = 11.68*8.85e-14;  % silicon permitivity F/cm
q = 1.602e-19;         % electron charge      
lam0 = c0/f0;

% update refractive index at fref
dlam = (c0/fref-c0/f0);
n0 = n0 + dlam*(n0-ng)/lam0;

% Soref & Bennett, 1987, unit: cm
A = lam0^2*3.64e-10;
B = lam0^2*3.51e-6;
C = 1/2*lam0^2*3.52e-6;  
D = 1/2*lam0^2*2.4e-6;
% 1/2 factor converting power coef. to field coef.


switch model_select
    case 1
        neff = n0 + n1*(-V);
        alpha = a0 + a1*(-V);
        
    case 2
        neff = n0 + n1*sqrt(Vbi-V);
        alpha = a0 + a1*sqrt(Vbi-V);
        
    case 3
        Vbi = 0.0259*log(Nd*Na/10^20);
        xp = sqrt(2*eps*Nd/(q*Na*(Nd+Na))*(Vbi-V));
        xn = sqrt(2*eps*Na/(q*Nd*(Nd+Na))*(Vbi-V));    
        % update n0, a0
        n0 = n0 - (A*Nd + B*Na^0.8)/2 * wf;
        a0 = a0 + 100*(C*Nd + D*Na)/2 * wf;
        neff = n0 + (xp*B*Na^0.8 + xn*A*Nd)/Lj * wf;
        alpha = a0 - 100*(xp*D*Na + xn*C*Nd)/Lj * wf;
        
    otherwise
        neff = n0;
        alpha = a0;
end
%#codegen
