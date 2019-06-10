close all; clearvars
dbstop if error

%ltx {\bf Input}\label{SBFEPoly2NInput}
probdef = @ProbDefLShapedPlate; %ltx handler of problem definition function \label{SBFEPoly2NProbDef}
%ltx mesh generation options
mesher = 1; %ltx = 1: DistMesh ; = 2 PolyMesher;\\ \hspace*{\fill } = 3: input triangle mesh; = 4: input ploygon mesh \label{SBFEPoly2NProbDef_mesher}
para.h0 = 0.2; %ltx element size when using DistMesh\label{SBFEPoly2NInput_h0}
para.nPolygon = 128; %ltx number of polygons when using PolyMesher\label{SBFEPoly2NInput_nPoly}
para.nDiv = 1; %ltx number of elements to divide one edge into

%ltx {\bf Mesh generation}\label{SBFEPoly2NMesh}
[ coord, sdConn, sdSC ] = createSBFEMesh(probdef, mesher, para);

figure; 
opt = struct('sdSC', sdSC, 'fill', [0.9 0.9 0.9]);
% % opt=struct('LineSpec','-b', 'sdSC', sdSC, 'LabelSC', 14, ...
% %     'fill', [0.9 0.9 0.9], 'PlotNode', 1, 'LabelNode', 12);
PlotSBFEMesh(coord, sdConn, opt); 
title('MESH');
% % %ltx eliminate round-off errors of mesher at boundary\label{SBFEPoly2NCoord}
% % %ltx it is only activated for the examples of rectangular domains
% % coord =  round(coord*1e5)/1e5;

% analysis
[U, sdSln] = SBFEPoly2NSolver(probdef, coord, sdConn, sdSC);