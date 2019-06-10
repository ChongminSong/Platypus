function [K, d, v11, M] = SElementSlnEigenMethod(E0, E1, E2, M0)
%Stiffness matrix of an S-Element
%
% Inputs:
%     E0, E1, E2, M0  - coefficient matrices of S-Element
%
% Outputs:
%     K      - stiffness matrix of S-Element
%     d      - eigenvalues 
%     v11    - upper half of eigenvectors (displacement modes)
%     M      - mass matrix of S-Element

nd = size(E0,1); % number of DOFs of boundary nodes
id = 1:nd;      % index of selected eigenvalues (working variable)

%ltx Preconditioning
Pf = 1./sqrt(abs(diag(E0))); P  = diag(Pf); %ltx Eq.~\eqref{eq:Precond-factor}
E0 = P*E0*P; E1 = P*E1*P; E2 = P*E2*P; %ltx Eq.~\eqref{eq:Precond-coeff}
%ltx Construct $[Z_p]$, see Eq.~\eqref{eq:Precond-Zp}
m = E0\[-E1' eye(nd)];
Zp = [m; E2+E1*m(:,id) -m(:,id)'];
%ltx eigenvalues and eigenvectors of $[\hat{Z}_p]$, see Eq.~\eqref{eq:Precond-eign}
[v, d] = eig(Zp); %ltx $\hat{v}$ : eigenvector matrix
d = diag(d); %ltx eigenvalues stored as a column vector
%ltx index for sorting eignvalues in descending order of real part.
[~, idx] = sort(real(d),'descend');
%ltx  select eigenvalues and eigenvectors for solution of bounded domain
d = d(idx(id)); %ltx select the first half of sorted eigenvalues
v = v(:, idx(id)); %ltx select the corresponding eigenvectors
v = diag([Pf;1./Pf])*v; %ltx Eq.~\eqref{eq:Precond-eign-vector}
%ltx modes of translational rigid body motion, see Item~\vref{enu:Zero-eigenvalues}
d(end-1:end) = 0; %ltx set last two eigenvalues to zero
v(:,end-1:end) = 0; %ltx set last two eigenvectors to zero
v(1:2:nd-1,end-1) = 1; %ltx set $u_{x}=1$ in $\{\phi^{\textrm{(u)}}_{n-1}\}$
v(2:2:nd,end) = 1; %ltx set $u_{y}=1$ in $\{\phi^{\textrm{(u)}}_{n}\}$
%ltx normalisation of eigenvectors, see Eq.~\eqref{eq:Precond-eign-vector-normalised}
v = bsxfun( @rdivide, v, max(abs(v(1:nd-2,:))) );
v11 = v(id, :); v11inv = inv(v11);
K  = real(v(nd+id, :)*v11inv); %ltx stiffness matrix, see Eq.~\eqref{eq:EigenMethod-stiffnessMatrix}

%ltx mass matrix; \label{Code_EigenMethod_mass}
if nargout > 3
    M0 = v11'*M0*v11; %ltx Eq.~\eqref{eq:2Dm0}
    %cbgnltx
    % \texttt{am} is a square matrix with all columns being the vector of
    % eigenvalues. The entry $(i,\ j)$ equals $\lambda_i$.
    % The entry $(i,\ j)$ of \texttt{am'} equals $\lambda_j$.
    %cendltx
    am = d(:,ones(1,nd));
    M0 = M0./(2+am+am'); %ltx the  entry $(i,\ j)$ of \texttt{am+am'} equals to $\lambda_i+\lambda_j$
    M = real(v11inv'*M0*v11inv); %ltx Eq.~\eqref{eq:2DM}
end
end