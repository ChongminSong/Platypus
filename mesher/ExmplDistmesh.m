clearvars; close all
dbstop if error

% probdef = @ProbDefPureBendingBeam;
% probdef = @ProbDefCantileverBeam; 
% probdef = @ProbDefPlateWithHole;
% probdef = @ProbDefLShapedPlate;
probdef = @ProbDefSquarePatch; %ltx specify problem definition function\label{examplesMeshProbDef}
mesher = 0; %ltx = 0: input triangle mesh; = 1 PolyMesher; = 2 distmesh2d
para.nPolygon = 128; %ltx number of polygons when using PolyMesher
para.h0 = 0.1; %ltx element size when using distmesh2d
[ coord, sdConn, sdSC ] = createSBFEMesh(probdef, mesher, para);
coord =  round(coord*1e6)/1e6; %ltx eliminate error on boundary

%ltx plot mesh
figure
PlotSBFEMesh(coord, sdConn);
title('MESH');

%ltx {\bf Materials}: elascity matrix and mass density
mat = probdef('MAT');

%ltx {\bf Boundary conditions}
F = zeros(2*size(coord,1),1);
glbMesh = struct('coord',coord , 'sdConn',{sdConn});
outArg = probdef('BCond', glbMesh);
BC_Disp = outArg{1}; %ltx displacement constrains
edge    = outArg{2}; %ltx nodes of edges with prescribed surface traction
trac    = outArg{3}; %ltx prescribed surface traction
if ~isempty(edge)
    F = addSurfTraction(coord, edge, trac, F); %ltx add equivalent nodal force
end

%ltx {\bf Static solution}
%ltx solution of subdomains and assemblage of global stiffness and mass matrices
[sdSln, K, M] = SBFEMAssemblage(coord, sdConn, sdSC, mat);
%ltx solution of nodal displacements and reactions
[U, ReactFrc] = SolverStatics(K, BC_Disp, F);
%ltx plot deformed mesh
figure
opt = struct('MagnFct', 0.1, 'Undeformed',':b');
PlotDeformedMesh(U, coord, sdConn, opt);
title('DEFORMED MESH');

output = probdef('Output', glbMesh);
disp('   DOF          Displacement');
fprintf('%6d  %25.15e\n',[output.DispDOF U(output.DispDOF)]');

ex = probdef('EXACT',coord);
if ~isempty(ex)
    Uex = reshape(ex(:,1:2)',[],1);
    err = Uex-U;
    eNorm = sqrt((err'*err)/(Uex'*Uex));
    disp([' Displacement error norm = ', num2str(eNorm)]);
end

% % strain modes of subdomains
% sdStrnMode =  SubdomainStrainMode2NodeEle( sdSln );
% % integration constants
% sdIntgConst = SubdomainIntgConst( U, sdSln );
% 
% isd = 1;
% xi = 1; % radial coordinate
% % displacements and strains at specified raidal coordinate
% [nodexy, dsp, strnNode, GPxy, strnEle] = ...
%     SubdomainInDispStrain(xi, sdSln{isd},  ...
%     sdStrnMode{isd}, sdIntgConst{isd});
% 
% exNode = probdef('EXACT',nodexy);
% exGP = probdef('EXACT',GPxy);
% 
% [exNode(:,3) (mat.D*strnNode)']
% [exGP(:,3) (mat.D*strnEle)']
% 
