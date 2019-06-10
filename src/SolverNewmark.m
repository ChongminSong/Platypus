function [t, dsp, vel, accl] = SolverNewmark(TIMEPara, ...
                              forceHistory, K, M, BC_Accl, F)
%Newmark time integrator 
%
%Inputs:
%   TIMEPara     : array containing parameters for Newmark method
%                  TIMEPara(1): number of time steps
%                  TIMEPara(2): size of time step 
%                  TIMEPara(3): gamma, optional (default is 0.5) 
%                  TIMEPara(4): beta, optional (default is 0.25)
%   forceHistory : external action at a given time is multiplied   
%                  with a factor obtain by interpolating the array 
%                  [0 f0; t1 f1; t2 f2, t3 f3, ...] 
%   K            : stiffness matrix
%   M            : mass matrix
%   BC_Accl      : prescribed nodal accelerations in rows of
%                  [Node, Dir. (=1 for x; =2 for y), Acceleration]
%   F            : load vector 
%
%Outputs:
%   t            : time 
%   dsp          : response of nodal displacements
%   vel          : response of nodal velocity
%   accl         : response of nodal acceleration

ns = TIMEPara(1); %ltx number of time steps
dt = TIMEPara(2); %ltx size of time step
if length(TIMEPara)>2 %ltx parameters $\gamma$ and $\beta$ 
    gamma = TIMEPara(3);
    beta =  TIMEPara(4);
else %ltx default values 
    gamma = 0.5;
    beta =  0.25;
end

ndn = 2; %ltx 2 DOFs per node
NDof = size(K,1); %ltx number of DOFs

%ltx displacement boundary condition
CDofs = []; %ltx initilising constrained DOFs
if ~isempty(BC_Accl)
    %ltx constrained DOFs with prescribed accelerations
    CDofs = (BC_Accl(:,1)-1)*ndn + BC_Accl(:,2);
end
FDofs   = (1:NDof); %ltx initilising unconstrained (free) DOFs
FDofs(CDofs) = []; %ltx remove constrained DOFs

t = (0:ns)*dt; %ltx time at all time station $0$, $\Delta_t$, $2\Delta_t$, $\ldots$
%ltx load factor at all time station by linear inter- or extra-polation
ft = interp1(forceHistory(:,1), forceHistory(:,2), t,...
             'linear', 'extrap' );

%ltx initialising variables storing response history for output        
dsp = zeros(NDof,ns+1); %ltx displacements
vel = dsp;  %ltx velocities
accl = dsp; %ltx accelerations

%ltx initialising variables
dn = zeros(NDof,1); %ltx displacements
vn = dn; %ltx velocities
an = dn; %ltx accelerations

%ltx initial time step (t = 0)
if abs(ft(1)) > 1.d-20
    %ltx balance initial forces with inertial forces
    f = ft(1)*F; %ltx external nodal forces at $n=0$
    an(CDofs) = ft(1)*BC_Accl(:,3); %ltx prescribed accelerations  
    if ~isempty(CDofs)
        %ltx enforcing prescribed accelerations, Eq.~\eqref{eq:Newmark-step0}
        f = f - M(:,CDofs)*an(CDofs); 
    end
    %ltx solution of accelerations at unconstrained DOFs, Eq.~\eqref{eq:Newmark-step0}
    an(FDofs) = (M(FDofs,FDofs)\f(FDofs));
end

%ltx wroking variables, Eqs.~\eqref{eq:Newmark-corrector} and \eqref{eq:Newmark-predictor} 
fg1 = gamma*dt;
fg2 = (1.d0-gamma)*dt;
fb1 = beta*dt*dt;
fb2 = (0.5d0-beta)*dt*dt;

%ltx dynamic-stiffness matrix of unconstrained DOFs 
%ltx $\beta\Delta_{t}^{2}[K_{11}]+[M_{11}]$, Eq~\eqref{eq:Newmark-stepn} 
dstf = M(FDofs,FDofs) + fb1*K(FDofs,FDofs);
%ltx Cholesky factorization of dynamic stiffness matrix
[Ldstf, ~, s] = chol(dstf,'lower','vector');

tmp = zeros(length(FDofs),1);
for it = 2:ns+1
    dn = dn + dt*vn + fb2*an; %ltx displacement predictor, Eq.~\eqref{eq:Newmark-predictor-disp}
    vn = vn + fg2*an;         %ltx velocity predictor, Eq.~\eqref{eq:Newmark-predictor-disp}
    f = ft(it)*F; %ltx external nodal forces
    if ~isempty(CDofs)
        an(CDofs) = ft(it)*BC_Accl(:,3); %ltx prescribed accelerations  
        %ltx correctors for constrained DOFs, Eq.~\eqref{eq:Newmark-predictor}
        vn(CDofs) = vn(CDofs) + fg1*an(CDofs); %ltx velocities
        dn(CDofs) = dn(CDofs) + fb1*an(CDofs); %ltx displacements
        %ltx add $-[M_{12}]\{\ddot{u}_{2}\}_{n} $ to right-hand side, Eq.~\eqref{eq:Newmark-stepn}
        f = f - M(:,CDofs)*an(CDofs);
    end
    %ltx add $-[K_{11}]\{\tilde{u}_{1}\}_{n}-[K_{12}]\{u_{2}\}_{n} $ to right-hand side, Eq.~\eqref{eq:Newmark-stepn}
    f = f - K*dn;
    f = f(FDofs); %ltx right-hand side of Eq.~\eqref{eq:Newmark-stepn}
    
    %ltx solve Eq.~\eqref{eq:Newmark-stepn} for unknown accelerations
    tmp(s) = Ldstf'\(Ldstf\(f(s)));
    an(FDofs) = tmp; %ltx update acceleration vector containing all DOFs
    %ltx correctors, Eq.~\eqref{eq:Newmark-corrector}
    vn(FDofs) = vn(FDofs) + fg1*an(FDofs);  %ltx displacements 
    dn(FDofs) = dn(FDofs) + fb1*an(FDofs);  %ltx velocities
    
    %ltx store responses for output
    dsp(:,it) = dn;
    vel(:,it) = vn;
    accl(:,it) = an;
    
end
