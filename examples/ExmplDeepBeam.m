%cbgnltx
%\leftadjust {\bf Matlab script: Cantilever beam under bending}  
%cendltx
%ltx {\bf Mesh}
%ltx nodal coordinates. One node per row \texttt{[x y]}\label{Exmpl_rectangle_polygon_node}
coord = [    0      1;    0    0.5;  0        0;
          0.68      1; 0.68   0.63;  0.49  0.48;
           0.5      0; 1.35      1;  2        1;
             1   0.45; 1.35   0.62;  1        0;
           1.5   0.47;    2    0.5;  1.5      0;
             2      0  ];
%ltx nodes of a ploygon. The sequence follows counter-clock direction.\label{Exmpl_rectangle_polygon}          
polygon = { [ 3     7     6     2];
            [15    16    14    13];
            [ 2     6     5     4     1];
            [12    15    13    11    10];
            [ 7    12    10     5     6];
            [ 4     5    10    11     8];
            [ 8    11    13    14     9] };
%cbgnltx
% Input S-element connectivity as a cell array (One S-element per cell). 
% In a cell, the connectivity of line elements is given by
% one element per row \texttt{[Node-1 Node-2]}.
%cendltx
nsd = length(polygon); %ltx number of S-elements
sdConn = cell(nsd,1);  %ltx initialising connectivity
sdSC = zeros(nsd,2);   %ltx scaling centre
for isub =1:nsd
    %ltx build connectivity \label{Exmpl_Rectangle_conn}
    sdConn{isub}=[polygon{isub}; ...
                  polygon{isub}(2:end) polygon{isub}(1)]';
    %ltx scaling centre at centroid of nodes (averages of nodal coorindates)\label{Exmpl_Rectangle_sc}
    sdSC(isub,:) = mean(coord(polygon{isub},:));
end

%ltx {\bf Materials}: elascity matrix for plane stress condition
mat.D = IsoElasMtrx(10E9, 0.2); %ltx  ($E$ in Pa, $\nu$) 
mat.den = 2000; %ltx mass density in $\unit{kg/m^{3}}$

%ltx {\bf Boundary conditions} 
%ltx displacement constrains (or prescribed acceleration in a response 
%ltx  history analysis). One constrain per row: \texttt{[Node Dir]}
BC_Disp = [1 1 0; 1 2 0; 2 1 0; 2 2 0; 3 1 0; 3 2 0];

%ltx assemblage nodal forces 
ndn = 2; %ltx 2 DOFs per node
NDof = ndn*size(coord,1); %ltx number of DOFs
F = zeros(NDof,1); %ltx initialising right-hand side of equation $[K]{u} = {F}$
edge = [16 14; 14 9]; %ltx edges subject to tractions, one row per edge
trac = [6E5 0 0 0; 0  0 -6E5 0]'; %ltx tractions, one column per edge
F = addSurfTraction(coord, edge, trac, F);
%ltx plot mesh
figure
opt=struct('sdSC',sdSC, 'LabelSC', 14,'LineSpec','-k', ...
           'PlotNode',1, 'LabelNode', 14,...
           'BC_Disp',BC_Disp);
PlotSBFEMesh(coord, sdConn, opt);
title('MESH');