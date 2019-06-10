%InfPlateWithHole
close all;clearvars;dbstop error;

% handler of problem definition function
probdef=@ProbDefInfPlateWithHole;

%% Mesh generation
mesher=2;
para.h0=0.05; % element size when using Distmesh
para.nPolygon=512*2;
[coord, sdConn, sdSC]=createSBFEMesh(probdef, mesher, para);

opt=struct('sdSC', sdSC, 'fill', [0.9 0.9 0.9]);
PlotSBFEMesh(coord, sdConn, opt);
title('MESH');

%% static analysis 
slnOpt.type='STATICS';
[U, sdSln]= SBFEPoly2NSolver(probdef,coord,sdConn,sdSC,slnOpt);

%% post processing
% compute stress at nodes and scaling centre
[SCStrs, NodeStrs]= SBFEPolyStress(probdef, U, sdSln, sdConn, sdSC); 
% plot stress contour
Poly2triMeshContour(coord, sdSC, SCStrs, NodeStrs, sdConn);
 