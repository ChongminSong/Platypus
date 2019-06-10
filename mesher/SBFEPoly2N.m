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

%ltx {\bf Material}: elascity matrix and mass density\label{SBFEPoly2NMaterial}
mat = probdef('MAT'); %ltx get material constants

%ltx {\bf Boundary conditions} \label{SBFEPoly2NBC}
%ltx build input argument of problem definition function
% % glbMesh = struct('coord',coord , 'sdConn',{sdConn});
% % outArg = probdef('BCond', glbMesh); %ltx get boundary conditions
outArg = probdef('BCond', {coord, sdConn}); %ltx get boundary conditions
%ltx interpret output argument of problem definition function
BC_Disp = outArg{1}; %ltx displacement constrains
F       = outArg{2}; %ltx nodal forces

%ltx {\bf Solution} \label{SBFEPoly2NSolution}
%ltx subdomain solution and global stiffness matrix
[sdSln, K, M] = SBFEMAssemblage(coord, sdConn, sdSC, mat);
%ltx nodal displacements and reactions
[U, ReactFrc] = SolverStatics(K, BC_Disp, F);

%ltx {\bf Post-processing} \label{SBFEPoly2NPostProcessing}
%ltx plot deformed mesh
figure
opt = struct('MagnFct', 0.1, 'Undeformed',':b');
PlotDeformedMesh(U, coord, sdConn, opt);
title('DEFORMED MESH');
%ltx output displacements at selected DOFs
output = probdef('Output', {coord, sdConn});
if ~isempty(output)
    disp('   DOF          Displacement');
    fprintf('%6d  %25.15e\n',[output.DispDOF U(output.DispDOF)]');
end
%ltx compute error norm
Uex = probdef('EXACT',coord);
if ~isempty(Uex)
    %ltx reshape the exact solution to a column vector
    Uex = reshape(Uex(:,1:2)',[],1);
    err = Uex-U; %ltx displacement error
    eNorm = norm(err)/norm(Uex);
    disp([' Displacement error norm = ', num2str(eNorm)]);
end