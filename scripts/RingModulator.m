%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sen Lin <senlin@berkeley.edu>
% Integrated Systems Group, EECS, UC Berkeley
% 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ring modulator class

classdef RingModulator
    properties
        % model select for phase shifter
        % 1--linear, 2--sqrt, 3--diode
        model = 1;      % default linear
        L = 30e-6;      % ring perimeter
        t_in = 0.98;    % input port coupling
        t_drop = 1;     % drop port coupling
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
        function ring = RingModulator(model, L, t_in, t_drop, ...
                n0, ng, a0, Ne, Nh, Lj, wf, a1, n1, Vbi)
            if nargin == 14
                ring.model = model;
                ring.L = L;
                ring.t_in = t_in;
                ring.t_drop = t_drop;
                ring.n0 = n0;
                ring.ng = ng;
                ring.a0 = a0;
                ring.Ne = Ne;
                ring.Nh = Nh;
                ring.Lj = Lj;
                ring.wf = wf;
                ring.a1 = a1;
                ring.n1 = n1;
                ring.Vbi = Vbi;             
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
        
        function [Et, Ed, Pt, Pd] = tf(self, lambda, V)
            % transfer function derived based on T-matrix method

            [neff, alpha] = self.material(lambda, V);
            
            % phase shift
            theta = 2*pi*self.L*neff./lambda;
            % round-trip absorption
            a_rt = exp(-self.L*alpha);
            
            % half ring phase shift and absorption
            a_h = sqrt(a_rt);
            theta_h = theta/2;
            
            % loseless coupler;
            k1 = 1i * sqrt(1-self.t_in^2);
            k2 = 1i * sqrt(1-self.t_drop^2);
            
            % drop port E-field transfer function
            Ed = - conj(k1)*k2*a_h.*exp(1i*theta_h)./...
                (1-conj(self.t_in).*conj(self.t_drop).*a_rt.*exp(1i*theta));  
            
            % output port E-field transfer function
            Et = (self.t_in - conj(self.t_drop).*a_rt.*exp(1i*theta))./...
                (1-conj(self.t_in).*conj(self.t_drop).*a_rt.*exp(1i*theta));
            
            % output port (Pt) and drop port (Pd) power
            Pt = abs(Et).^2;
            Pd = abs(Ed).^2;
            
        end
        

        function [laser_opt, OMA_max, P1, P0] = modulation(self, lambda, V0, V1)
            % find the optimal laser wavelength (relative to resonance) to
            % maximize OMA
            [~,~,Pt0,~] = self.tf(lambda, V0);
            [~,~,Pt1,~] = self.tf(lambda, V1);
            OMA = Pt1 - Pt0;
            [OMA_max,ind] = max(OMA);
            P1 = Pt1(ind);
            P0 = Pt0(ind);
            laser_opt = lambda(ind);
            
        end
        
        function [state, t_drop_critical] = coupling_state(self,lambda0, V) 
            % find the coupling state of the ring, assuming absorption
            % changes with voltage bias
            % find the nominal drop port coupling to
            % get closest to critical coupling
            [~, alpha] = self.material(lambda0, V);
            if self.t_in^2 > self.t_drop^2*exp(-self.L*alpha*2)
                state = 'under coupled';
                t_drop_critical = 1;
            elseif self.t_in^2 == self.t_drop^2*exp(-self.L*alpha*2)
                state = 'critically coupled';
                t_drop_critical = self.t_drop;
            else
                state = 'over coupled';
                t_drop_critical = sqrt((self.t_in^2)/exp(-self.L*alpha*2));
            end
        end
            
    end 
end