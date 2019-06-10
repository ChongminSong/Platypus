function [x] = ProbDefDeepBeam(Demand,Arg)
%Problem definition function of a deep beam
%
%inputs
%   Demand    - option to access data in local functions
%   Arg       - argument to pass to local functions
%output
%   x         - cell array of output.

%ltx \texttt{L}: Length; \texttt{H}: Height;
%ltx \texttt{E}: Young's modulus; \texttt{pos}: Poisson's ratio; \texttt{den}: mass density
para = struct('L',2, 'H',1, 'M', 100, ...
    'E', 10E6, 'pos',0.2, 'den',2);
%ltx bounding box [x1,x2,y1,y2]
BdBox = [0 para.L -para.H/2 para.H/2];
switch(Demand) %ltx access local functions according to \texttt{Demand}\label{Code_ProbDefPureBendingBeam_case}
    case('para');    x = para; %ltx parameters of problem definitions
        
    case('TriMesh'); x = TriMesh; %ltx direct input triangular mesh
    case('RectMesh'); x = RectMesh(Arg, para);
        
    %ltx for use with polygon mesh generator PolyMesher \label{Code_ProbDefPureBendingBeam_PolyMesher}
    case('Dist');     x = DistFnc(Arg,BdBox);
    case('BC');       x = {[],[]};
        
    %ltx for use with triangular mesh generator DistMesh \label{Code_ProbDefPureBendingBeam_DistMesh}
    case('Dist_DistMesh');  x = DistFnc(Arg,BdBox);
        x = x(:,end); %ltx only use the last column
    case('fh');     x = huniform(Arg);
    case('pfix');   x = [BdBox([1 3]); BdBox([1 4]); ...
            BdBox([2 3]); BdBox([2 4])];
        
    case('BdBox');  x = BdBox; %ltx for use with PolyMesher and DistMesh
   
    case('MAT');     x = struct('D',IsoElasMtrx(para.E, ...
            para.pos), 'den', para.den); %ltx material constants
    
    case('BCond');   x = BoundaryCond(Arg, para, BdBox);
    case('BCondModal');   x = BoundaryCondModal(Arg, para, BdBox);
    case('BCondTime');   x = BoundaryCondTime(Arg, para, BdBox);
    
    case('Output');  x = OutputRequests(Arg, BdBox);
    case('OutputTime');  x = OutputRequests(Arg, BdBox);
    
    case('EXACT');   x = ExactSln(Arg, para); %ltx exact solution
    otherwise
        warning('Unexpected keyword in ProbDefDeepBeam.')
end
end

%ltx \textbf{Signed distance function}\label{Code_ProbDefPureBendingBeam_DistFnc}
function Dist = DistFnc(P,BdBox)
Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
end

%ltx \textbf{Boundary conditions}  \label{Code_ProbDefPureBendingBeam_BC}
function [x] = BoundaryCond(Arg, beam, BdBox)
coord = Arg{1}; %ltx nodal coordinates
sdConn = Arg{2}; %ltx element connectivity

eps = 1.d-4*sqrt(beam.L*beam.H/size(coord,1));
%ltx \emph{Displacement constrains} (prescribed acceleration in a response history 
%ltx analysis). One constrain per row: \texttt{[Node Dir]}
lNodes = find(abs(coord(:,1)-BdBox(1))<eps); %ltx nodes at left edge
llNode = find(abs(coord(lNodes,1)-BdBox(1))<eps & ...
    abs(coord(lNodes,2)-BdBox(3))<eps); %ltx lower left corner
n1 = length(lNodes);
ex = ExactSln(coord(lNodes,:), beam); %ltx exact solution at left side
BC_Disp = [ [lNodes; lNodes(llNode)], [ones(n1,1); 2], ...
    [ex(:,1); ex(llNode,2)]]; %ltx displacement boundary condition

%ltx \emph{Surface traction at right side of beam}
%ltx edges of polygon S-elements
MeshEdges = meshConnectivity( sdConn );
%ltx centres of edges of S-elements
centres = (coord(MeshEdges(:,1),:)+coord(MeshEdges(:,2),:))/2;
%ltx find edges at the right side of the beam
rEdges = MeshEdges(abs(centres(:,1)-BdBox(2))<eps,:);
%ltx surface traction at nodes
ex = ExactSln(coord([rEdges(:,1); rEdges(:,2)],:), beam);
n2  = length(rEdges);
trac = [ex(1:n2,3)'; zeros(1,n2); ex(n2+1:end,3)'; zeros(1,n2)]; 
%ltx nodal forces
F = zeros(2*size(coord,1),1); %ltx initialising force vector
%ltx add nodal forces equivalent to surface traction
F = addSurfTraction(coord, rEdges, trac, F);

%ltx output
x = {BC_Disp, F};
end

%ltx \textbf{Boundary conditions for modal analysis}  \label{Code_ProbDefPureBendingBeam_BCModal}
function [x] = BoundaryCondModal(Arg, beam, BdBox)
coord = Arg{1}; %ltx nodal coordinates

eps = 1.d-4*sqrt(beam.L*beam.H/size(coord,1));
lNodes = find(abs(coord(:,1)-BdBox(1))<eps); %ltx nodes at left edge
n1 = length(lNodes);
BC_Disp = [ [lNodes; lNodes], [ones(n1,1); 2*ones(n1,1)], ...
    zeros(2*n1,1)]; %ltx displacement boundary condition

%ltx output
x = {BC_Disp};
end

%ltx \textbf{Boundary conditions for response history analysis}  \label{Code_ProbDefPureBendingBeam_BCTime}
function [x] = BoundaryCondTime(Arg, beam, BdBox)
coord = Arg{1}; %ltx nodal coordinates
sdConn = Arg{2}; %ltx element connectivity

eps = 1.d-4*sqrt(beam.L*beam.H/size(coord,1));
lNodes = find(abs(coord(:,1)-BdBox(1))<eps); %ltx nodes at left edge
n1 = length(lNodes);
BC_Disp = [ [lNodes; lNodes], [ones(n1,1); 2*ones(n1,1)], ...
    zeros(2*n1,1)]; %ltx displacement boundary condition

%ltx surface traction at the top of the beam
%ltx edges of polygon S-elements
MeshEdges = meshConnectivity( sdConn );
%ltx centres of edges of S-elements
centres = (coord(MeshEdges(:,1),:)+coord(MeshEdges(:,2),:))/2;
%ltx find edges at the top of the beam
tEdges = MeshEdges(abs(centres(:,1)-BdBox(2))<eps,:);
F = zeros(2*size(coord,1),1); %ltx initialising force vector
trac = [ 0; 100; 0; 100];  %ltx uniform surface traction
F = addSurfTraction(coord, tEdges, trac, F); %ltx equivalent nodal forces 
forceHistory = [0 0; 0.01 1; 0.03 -1; 0.04 0; 100 0];

%ltx output
x = {BC_Disp, F, forceHistory};
end

%ltx {\bf Output Requests}
function [outputs] =OutputRequests(Arg, BdBox)
coord = Arg{1}; %nodal coordinates
%ltx find the node at upper-right corner
inode = find( abs(coord(:,1)-BdBox(2)) <1.d-5 & ...
    abs(coord(:,2)-BdBox(4)) <1.d-5 );
outputs = struct('DispDOF', 2*inode); %ltx request vertical displacement
end

%ltx {\bf Exact solution} (see Eq.~\eqref{eq:PureBendingBeamExact})
function [x] = ExactSln(xy, beam)

x = xy(:,1);  y = xy(:,2); %ltx coordinates
M = beam.M; %ltx bending moment
I = beam.H^3/12; %ltx moment of inertia
E = beam.E; pos = beam.pos;
ux = -M/(E*I).*x.*y; %ltx displacements
uy = 0.5*M/(E*I).*( x.^2 + pos*(y.^2-(beam.H/2)^2) );
% % uy = 0.5*M/(E*I).*( x.^2 + pos*y.^2 );
sx = -M*y/I; sy = zeros(length(x),1); sxy = sy; %ltx stresses
x  = [ux, uy, sx, sy, sxy]; %ltx output
end

%ltx {\bf A pre-stored triangular mesh} \label{Code_ProbDefPureBendingBeam_trimesh}
function [x] = TriMesh()
%ltx a triangular mesh

%ltx nodal coordinates
p = [ 0.00  -0.16;   0.00  0.16;  0.00 -0.50;  0.00  0.50;
      0.41  -0.50;   0.41  0.50;  0.51 -0.00;  0.80 -0.50;
      0.80   0.50;   1.00  0.00;  1.20  0.50;  1.20 -0.50;
      1.49   0.00;   1.59  0.50;  1.59 -0.50;  2.00 -0.50;
      2.00   0.50;   2.00  0.16;  2.00 -0.16 ];

%ltx triangular elements
t = [5   1   3;  15  16  19;   7   2   1;   1    5   7;
     5   8   7;  18  17  14;   4   2   6;   2    7   6;
    10   8  12;  10   7   8;  12  15  13;  13   10  12;
    13  15  19;  19  18  13;  13  18  14;   9    6   7;
     7  10   9;  11   9  10;  11  13  14;  10   13  11];
%ltx output
x = {p, t};
end

%ltx {\bf Rectangular mesh} \label{Code_ProbDefPureBendingBeam_rectmesh}
function [x] = RectMesh(para, beam)
%ltx a rectangular mesh

h0 = para.h0; %ltx desired size of square S-elements
L = beam.L; H = beam.H; %ltx dimensions
nx = ceil(L/h0); ny = ceil(H/h0); %ltx number of S-elements
dx = L/nx; dy = H/ny; %ltx actual size of S-elements
[x, y] = meshgrid((0:nx)*dx, (0:ny)*dy-H/2); %ltx grid of nodes
coord = [x(:), y(:)]; %ltx nodal coordinates

%ltx construct  S-elements
[ic, ir] = meshgrid(0:nx-1, 0:ny-1); %ltx grid of S-elements
i0 = (ny+1)*(ic(:))+ir(:)+1; %ltx node at low-left corner of an S-element
sdConn = cell(nx*ny,1);  %ltx initialising connectivity
for ii = 1:nx*ny %ltx loop over S-elements
    i1 = i0(ii);  i2 = i1+ny+1; % lower left and right nodes
    sdConn{ii}=[i1 i2; i2 i2+1; i2+1 i1+1; i1+1 i1];
end
%ltx scaling centre
[x, y] = meshgrid((0.5:nx)*dx, (0.5:ny)*dy-H/2); %ltx grid of points
sdSC = [x(:), y(:)]; %ltx coordinates of scaling centres

%ltx output
x = struct('coord',coord , 'sdConn',{sdConn}, 'sdSC',sdSC );
end