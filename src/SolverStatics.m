function [d, F] = SolverStatics(K, BC_Disp, F)
%2D Static Analysis by the Scaled Boundary Finite Element Method
%
%Inputs:
%  K        - static stiffness matrix
%  BC_Disp  - prescribed displacements in rows of 
%             [Node, Direction (=1 for x; =2 for y), Displacement]
%  F        - external load vector 
%
%Outputs:
%  d        - nodal displacements
%  F        - external nodal forces including support reactions

ndn = 2; %ltx 2 DOFs per node
NDof = size(K,1);
d = zeros(NDof,1); %ltx Initialisation of nodal displacements

%ltx enforcing displacement boundary condition
%ltx initialisation of unconstrained (free) DOFs with unknown displacements
FDofs  = 1:NDof; 
if ~isempty(BC_Disp)
    %ltx constrained DOFs with prescribed displacements
    CDofs = (BC_Disp(:,1)-1)*ndn + BC_Disp(:,2);
    FDofs(CDofs) = []; %ltx remove constrained DOFs
    F = F - K(:,CDofs)*BC_Disp(:,3); %ltx Eq.~\eqref{eq:GlobalEquationSln-2}
    d(CDofs) = BC_Disp(:,3); %ltx store prescribed displacements
end

%ltx displacement of free DOFs, see Eq.~\eqref{eq:GlobalEquationSln-3}
d(FDofs) = K(FDofs,FDofs)\F(FDofs);
%ltx external forces, see Eq.~\eqref{eq:GlobalEquationSln-1}
F = K*d;

end