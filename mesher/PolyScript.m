clearvars; 
close all
dbstop if error

probdef = @PlateWithHoleDomain; %ltx specify function defining the model
% probdef = @CantileverDomain; %ltx specify function defining the model
% probdef = @CantileverPureBendingDomain;
%ltx creat polygon mesh
[coord, polygon] = PolyMesher(probdef,320,100); 
%convert polygons into subdomains of the SBFEM
[sdConn, sdSC] = polygonToSBFESubdomain(coord, polygon);

%ltx plot mesh
figure
opt=struct('LineSpec','-b', 'sdSC',sdSC, ...
     'PlotNode',0, 'LabelNode', 0,...
      'BC_Disp',[]); %ltx plotting options
%opt = [];
PlotSBFEMesh(coord, sdConn, opt);
title('MESH');

%ltx {\bf Materials}: elascity matrix for plane stress condition
%ltx \texttt{E}: Young's modulus; \texttt{p}: Poisson's ratio
ElasMtrx = @(E, p) E/(1-p^2)*[1 p 0;p 1 0;0 0 (1-p)/2];
outArg = probdef('MAT');
mat.D = ElasMtrx(outArg{1}, outArg{2}); %ltx  \texttt{E} in KPa
mat.den = outArg{3}; %ltx mass density in $\unit{Mg/m^{3}}$

%ltx {\bf Boundary conditions} 
%ltx initialising right-hand side of equation $[K]{u} = {F}$
F = zeros(2*size(coord,1),1); 

glbMesh = struct('coord',coord , 'sdConn',{sdConn});
outArg = probdef('myBC', glbMesh);
BC_Disp  =  outArg{1}; %ltx displacement constrains
edge = outArg{2}; %ltx nodes of edges with prescribed surface traction
trac    = outArg{3}; %ltx prescribed surface traction

if ~isempty(edge) 
    F = addSurfTraction(coord, edge, trac, F); %ltx add equivalent nodal force 
end
% static solution
SBFEM2DStaticAnalysisScript

ex = probdef('EXACT',coord);
Uex = reshape(ex(:,1:2)',[],1);
err = Uex-U;
sqrt(err'*err)/sqrt(Uex'*Uex)
[sqrt(numel(coord)) sqrt(size(coord,1)) sqrt(length(sdConn))]
sum(F)
inode = find(abs(coord(:,1)-max(coord(:,1)))<0.00001...
    &abs(coord(:,2)-min(coord(:,2)))<0.00001);
(Uex(2*inode)-U(2*inode))/Uex(2*inode)

% strain modes of subdomains
sdStrnMode =  SubdomainStrainMode2NodeEle( sdSln );

% integration constants 
sdIntgConst = SubdomainIntgConst( U, sdSln );

isd = 135;
xi = 1; % radial coordinate
% displacements and strains at specified raidal coordinate
[nodexy, dsp, strnNode, GPxy, strnEle] = ...
      SubdomainInDispStrain(xi, sdSln{isd},  ...
                            sdStrnMode{isd}, sdIntgConst{isd});

exNode = probdef('EXACT',nodexy);
exGP = probdef('EXACT',GPxy);

[exNode(:,3) (mat.D*strnNode)']
[exGP(:,3) (mat.D*strnEle)']
