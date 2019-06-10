%cbgnltx
%\leftadjust {\bf Matlab script: Uniaxial tension of a rectanglur plate 
% modelled with 3 subdomains}  
%cendltx
clearvars; close all; 
dbstop if error

%ltx {\bf Mesh}
%ltx nodal coordinates. One node per row \texttt{[x y]}
coord = [0 0; 0  1; 0 3; 1 0; 1 1; 2 0; 2 1; 2 3];
%cbgnltx
% Input subdomain connectivity as a cell array (One subdomain per cell). 
% In a cell, the connectivity of line elements is given by
% one element per row \texttt{[Node-1 Node-2]}.
%cendltx
sdConn = { [2 5; 5 7; 7 8; 8 3; 3 2]; %ltx subdomain 1
           [1 4; 4 5; 5 2; 2 1];      %ltx subdomain 2
           [4 6; 6 7; 7 5; 5 4]};     %ltx subdomain 3
%ltx coordinates of scaling centres of subdomains.
sdSC = [ 1 2; 0.5 0.5; 1.5 0.5]; %ltx one subdomain per row

%ltx {\bf Materials}
%ltx \texttt{E}: Young's modulus; \texttt{p}: Poisson's ratio
ElasMtrx = @(E, p) E/(1-p^2)*[1 p 0;p 1 0;0 0 (1-p)/2];
mat.D = ElasMtrx(10E9, 0.25); %ltx  \texttt{E} in Pa
mat.den = 2000; %ltx mass density in $\unit{kg/m^{3}}$\label{Exmpl_Uniaxial3SD_Sln_endMaterial}

%ltx {\bf Boundary conditions} \label{Exmpl_Uniaxial3SD_Sln_BC}
%ltx displacement constrains. One constrain per row: \texttt{[Node Dir Disp]}
BC_Disp = [1 2 0; 4 1 0; 4 2 0; 6 2 0];
%ltx nodal forces. One node per row: \texttt{[Node Dir F]} 
BC_Frc = [3  2  1E6; 8  2  1E6];  %ltx forces in N

%ltx{\bf Assemblage external forces} 
ndn = 2; %ltx 2 DOFs per node
NDof = ndn*size(coord,1); %ltx number of DOFs
F = zeros(NDof,1); %ltx initialising right-hand side of equation $[K]{u} = {F}$
F = AddNodalForces(BC_Frc, F); %ltx add prescribed nodal forces 

%ltx solution of subdomains and global stiffness and mass matrices
[sdSln, K, M] = SBFEMAssemblage(coord, sdConn, sdSC, mat);

%ltx {\bf Static solution}
dofs = (1:numel(coord))'; %ltx vector of DOFs for output 
SBFEM2DStaticsScript
