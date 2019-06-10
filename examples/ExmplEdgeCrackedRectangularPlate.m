%cbgnltx
%\leftadjust {\bf Matlab script: An edge-cracked rectangular body
% subject to tension}
%cendltx
clearvars; close all; dbstop if error

%ltx {\bf Mesh}
%ltx nodal coordinates in mm. One node per row \texttt{[x y]}
b = 0.1;
coord = b*[0 0; 0 -0.5; 0 -1; 0.5 -1; 1 -1; 1.5 -1; 2 -1;
    2 -0.5; 2 0; 2 0.5; 2 1; 1.5 1; 1 1; 0.5 1; 0 1; 0 0.5; 
    0 1E-14; 0 -2; 0 -3; 1 -3; 2 -3; 2 -2; 3 -3; 4 -3; 
    4 -2; 4 -1; 3 -1; 4 0; 4 1; 3 1; 4 2; 4 3; 3 3;
    2 3; 2 2; 1 3; 0 3; 0 2];
%cbgnltx
% Input S-element connectivity as a cell array (One S-element per cell).
% In a cell, the connectivity of line elements is given by
% one element per row \texttt{[Node-1 Node-2]}.
%cendltx
sdConn = { [1:16; 2:17]'; %ltx S-element 1
    [3 18; 18 19; 19 20; 20 21; 21 22; ...
    22 7;  4  3;  5  4;  6  5;  7  6];      %ltx S-element 2
    [21 23; 23 24; 24 25; 25 26; 26 27; ...
    27  7; 22 21;  7 22];  %ltx S-element 3
    [8  7;  9  8; 27 26;  7 27; 26 28; ...
    28 29; 29 30; 30 11; 10  9; 11 10];  %ltx S-element 4
    [30 29; 11 30; 29 31; 31 32; 32 33; ...
    33 34; 34 35; 35 11];  %ltx S-element 5
    [12 11; 13 12; 14 13; 15 14; 35 34; ...
    11 35; 34 36; 36 37; 37 38; 38 15]};  %ltx S-element 6
%ltx coordinates of scaling centres of S-elements, one S-element per row
sdSC = b*[ 1  0; 1 -2; 3 -2; 3 0; 3 2; 1 2]; 

%ltx {\bf Materials}
%ltx \texttt{E}: Young's modulus; \texttt{p}: Poisson's ratio
ElasMtrx = @(E, p) E/(1-p^2)*[1 p 0;p 1 0;0 0 (1-p)/2];
mat.D = ElasMtrx(10E9, 0.25); %ltx  \texttt{E} in Pa
mat.den = 2000; %ltx mass density in $\unit{kg/m^{3}}$

%ltx {\bf Boundary conditions} 
%ltx displacement constrains. One constrain per row: \texttt{[Node Dir Disp]}
BC_Disp = [19 1 0; 19 2 0; 20 1 0; 20 2 0; ...
      21 1 0; 21 2 0; 23 1 0; 23 2 0; 24 1 0; 24 2 0];
%ltx assemble load vector
ndn = 2; %ltx 2 DOFs per node
NDof = ndn*size(coord,1); %ltx number of DOFs
F = zeros(NDof,1); %ltx initialising right-hand side of equation $[K]\{u\} = \{F\}$
edge = [ 32 33; 33 34; 34 36; 36 37]; %ltx edges subject to traction
trac = [ 0 1E6 0 1E6]'; %ltx all edges have the same traction (in Pa),   
F = addSurfTraction(coord, edge, trac, F);

%ltx {\bf Plot mesh}
figure
opt=struct('sdSC',sdSC,'LineSpec','-k', ...
    'PlotNode',1, 'LabelNode', 1,...
    'BC_Disp',BC_Disp);
PlotSBFEMesh(coord, sdConn, opt);
title('MESH');