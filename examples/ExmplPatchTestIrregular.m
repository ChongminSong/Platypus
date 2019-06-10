%cbgnltx
%\leftadjust {\bf Patch test: A square modelled with 2 irregular polygon S-elements}
%cendltx
close all; clearvars; dbstop error;

%ltx {\bf Mesh}
%ltx nodal coordinates. One node per row \texttt{[x y]}
x1 = 0.5; y1 = 0.5; % Fig. b
% x1 = 0.05; y1 = 0.95; % Fig. c
coord = [ x1 y1; 0 0; 0.1  0; 1 0; 1 1; 0 1; 0 0.1];
%cbgnltx
% Input S-element connectivity as a cell array (One S-element per cell).
% In a cell, the connectivity of line elements is given by
% one element per row \texttt{[Node-1 Node-2]}.
%cendltx
sdConn = { [1 7; 7 2; 2 3; 3 1]; %ltx S-element 1
    [1 3; 3 4; 4 5; 5 6; 6 7; 7 1]};      %ltx S-element 2
%ltx coordinates of scaling centres of S-elements.\label{ExmplPatchTestIrregularSC}
if x1 > y1 %ltx extension of line $\overline{21}$ intersecting right edge of the square
    sdSC = [ x1/2 y1/2 ; (1+x1)/2 y1*(1+x1)/(2*x1)]; 
else %ltx extension of line $\overline{21}$ intersecting top edge of the square
    sdSC = [ x1/2 y1/2 ;  x1*(1+y1)/(2*y1) (1+y1)/2];
end
%ltx {\bf Materials}
mat.D = IsoElasMtrx(1, 0.25); %ltx elasticity matrix 
mat.den = 2; %ltx mass density 

%ltx {\bf Boundary conditions} \label{ExmplPatchTestIrregularBC}
%ltx displacement constrains. One constrain per row: \texttt{[Node Dir Disp]}
BC_Disp = [2 1 0; 2 2 0; 4 2 0];
%ltx assemble load vector
ndn = 2; %ltx 2 DOFs per node
NDof = ndn*size(coord,1); %ltx number of DOFs
F = zeros(NDof,1); %ltx initialising right-hand side of equation $[K]\{u\} = \{F\}$
% %ltx  horizontal tension 
% edge = [ 4 5; 6 7; 7 2]; 
% trac = [1 0 1 0; -1 0 -1 0; -1 0 -1 0]'; 
% %ltx vertical tension 
% edge = [2 3; 3 4; 5 6];
% trac = [0 -1 0 -1; 0 -1 0 -1; 0 1  0 1]'; 
%ltx pure shear
edge = [2 3; 3 4; 4 5; 5 6; 6 7;
        7 2]; %ltx edges subject to tractions, one row per edge
trac = [-1 0 -1 0; -1 0 -1 0; 0 1  0 1; 1 0 1 0;
        0 -1 0 -1; 0 -1 0 -1]'; %ltx tractions, one column per edge
F = addSurfTraction(coord, edge, trac, F);

%ltx {\bf Plot mesh}
figure
opt=struct('LineSpec','-k', 'sdSC',sdSC, ...
    'PlotNode',1, 'LabelNode', 1); %ltx plotting options
PlotSBFEMesh(coord, sdConn, opt);
title('MESH');

%ltx {\bf Static solution}
%ltx nodal displacements
SBFEM2DStaticAnalysisScript
disp('Nodal displacements')
for ii = 1:length(coord)
    fprintf('%5d %25.15e %25.15d\n',ii, d(2*ii-1:2*ii))
end

%ltx {\bf Stresses} \label{ExmplPatchTestIrregularStrs}
%ltx strain modes of S-elements
sdStrnMode =  SElementStrainMode2NodeEle( sdSln );
%ltx integration constants 
sdIntgConst = SElementIntgConst( d, sdSln );
%ltx displacements and strains at specified raidal coordinate
isd = 2; %ltx S-element number
xi = 1; %ltx radial coordinate
[nodexy, dsp, strnNode, GPxy, strnEle] = ...
      SElementInDispStrain(xi, sdSln{isd},  ...
                            sdStrnMode{isd}, sdIntgConst{isd});
disp('Stresses of Elements 1 and 2')
mat.D*strnEle(:,1:2)