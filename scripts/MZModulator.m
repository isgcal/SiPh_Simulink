%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sen Lin <senlin@berkeley.edu>
% Integrated Systems Group, EECS, UC Berkeley
% 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Mach Zenhder modulator class

classdef MZModulator
    properties
        % model select for phase shifter
        % 1--linear, 2--sqrt, 3--diode
        model = 1;      % default linear
        L = 30e-6;      % MZM arm length
        n0 = 2;         % effective index at f0
        ng = 3;         % group index
        a0 = 300;       % absorption coefficient (/m)
        Ne = 5e17;      % N-doping (/cm^3)
        Nh = 5e17;      % P-doping (/cm^3)
        Lj = 0.5e-4;    % junction geometry factor (cm)
        wf = 0.9;       % waveguide factor
        a1 = -0.3;      % alpha 1st order coef. for linear/sqrt model
        n1 = 8e-5;      % neff 1st order coef. for linear/sqrt model
        Vbi = 0.8;      % built-in voltage, only for sqrt model
        f0 = 230e12;    % measurement frequency for n0, a0
    end
    methods
        % constructor
        function mzm = MZModulator(model, L, n0, ng, a0, Ne, Nh, ...
                Lj, wf, a1, n1, Vbi)
            if nargin == 12
                mzm.model = model;
                mzm.L = L;
                mzm.n0 = n0;
                mzm.ng = ng;
                mzm.a0 = a0;
                mzm.Ne = Ne;
                mzm.Nh = Nh;
                mzm.Lj = Lj;
                mzm.wf = wf;
                mzm.a1 = a1;
                mzm.n1 = n1;
                mzm.Vbi = Vbi; 
            else
                disp('wrong number of inputs');
            end
        end
        
        function [neff, alpha] = material(self, lambda, V)
            % get the material properties for lambda and V
            c0 = 3e8;
            fref = c0./lambda;
            val= [self.model, self.n0, self.ng, self.a0, self.n1, ...
                self.a1,self.Ne,self.Nh,self.Vbi,self.Lj, self.f0, self.wf];
            len_f = length(fref);
            
            % material properties
            neff = zeros(1, len_f);
            alpha = zeros(1, len_f);
            
            for j = 1:len_f
                [neff(j), alpha(j)] = fmat(fref(j), V, val, c0);
            end        
        end
        
        function [Et, Ed, Pt, Pd] = tf(self, lambda, phase_os0, phase_os1, V0, V1)
            % transfer function derived based on T-matrix method
            % V0 and V1 are biases on each arm
            % phase_os is the offset phase on arm1
            [neff0, alpha0] = self.material(lambda, V0);
            [neff1, alpha1] = self.material(lambda, V1);
            % phase shift
            theta0 = 2*pi*self.L*neff0./lambda;
            theta1 = 2*pi*self.L*neff1./lambda;
            % single arm absorption
            a_arm0 = exp(-self.L*alpha0);
            a_arm1 = exp(-self.L*alpha1);
            
            % output 0 E-field transfer function
            Et = 0.5*exp(1i*(theta0+phase_os0)).*a_arm0 ...
            - 0.5 *exp(1i*(theta1+phase_os1)).*a_arm1; 
            Ed = 0.5*exp(1i*(theta0+phase_os0)).*a_arm0 ...
            + 0.5 *exp(1i*(theta1+phase_os1)).*a_arm1; 
            % output power at port0 and port1
            Pt = abs(Et).^2;
            Pd = abs(Ed).^2;

        end
        

        function [laser_opt, OMA_max, P1, P0] = modulation(self, lambda, phase_os0, phase_os1, V0, V1)
            % find the optimal laser wavelength to
            % maximize OMA
            [~,~,Pt0,~] = self.tf(lambda, phase_os0, phase_os1, V0, V1);
            [~,~,Pt1,~] = self.tf(lambda, phase_os0, phase_os1, V1, V0);
            OMA = Pt1 - Pt0;
            [OMA_max,ind] = max(OMA);
            P1 = Pt1(ind);
            P0 = Pt0(ind);
            laser_opt = lambda(ind);
            
        end
        
            
    end 
end