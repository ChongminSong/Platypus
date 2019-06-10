%cbgnltx
%\leftadjust  {\bf Matlab script: Uniaxial tension of a rectanglur plate 
% modelled with 3 subdomains}  
%cendltx
clearvars; close all; 

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
mat.den = 2000; %ltx mass density in $\unit{kg/m^{3}}$

%ltx {\bf Global stiffness matrix}
[sdSln, K, M] = SBFEMAssemblage(coord, sdConn, sdSC, mat);
