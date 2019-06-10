function [freq, modeShape] = SolverMode(NMode, K, M, BC_Disp)
%Compute natrural frequencies and modal shapes

ndn = 2; %ltx number of DOFs per node

NDof = size(K,1); %ltx number of DOFs
%ltx find unconstrained (free) DOFs from displacement boundary condition
%ltx initialisation of free DOFs with unknown displacements
FDofs = 1:NDof;
if ~isempty(BC_Disp)
    %ltx constrained DOFs with prescribed displacements
    CDofs = (BC_Disp(:,1)-1)*ndn +BC_Disp(:,2); 
    FDofs(CDofs) = []; %ltx remove constrained DOFs
end

%ltx solve eigenproblem considering unconstrained DOFs only, see Eq.~\eqref{eq:freeVibration-eigen}  
if NDof > 20*NMode && NDof > 500 %ltx compute the lowest \texttt{NMode} modes
    [modeShape, freq] = eigs(K(FDofs,FDofs),...
                               M(FDofs,FDofs), NMode, 'sm');
else %ltx compute all modes 
    [modeShape, freq] = eig(full(K(FDofs,FDofs)),...
                              full(M(FDofs,FDofs)));
end

%ltx sort eigefrequencies (square-root of eigenvalues) in ascending order
[freq, idx] = sort(sqrt(diag(freq)),'ascend'); %ltx Eq.~\eqref{eq:freeVibration-NatrualFreq}
%ltx arrange eigenvectors in the same order
mtmp = modeShape(:,idx(1:NMode)); 
%ltx initialise modal shapes to include constrained DOFs
modeShape = zeros(NDof,NMode);
modeShape(FDofs,:) = mtmp; %ltx fill in unconstrained DOFs 
freq = freq(1:NMode); %ltx output selected natural frequencies