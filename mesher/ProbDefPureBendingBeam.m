function [x] = ProbDefPureBendingBeam(Demand,Arg)
%Problem definition function of pure bending of a beam
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
    case('TriMesh'); x = TriMesh; %ltx direct input of a triangular mesh
    case('RectMesh'); x = RectMesh(Arg, para);
        %ltx mesh generators
        %ltx for use with PolyMesher \label{Code_ProbDefPureBendingBeam_PolyMesher}
    case('Dist');     x = DistFnc(Arg,BdBox);
    case('BC');       x = {[],[]};
        
        %ltx for use with DistMesh \label{Code_ProbDefPureBendingBeam_DistMesh}
    case('Dist_DistMesh');  x = DistFnc(Arg,BdBox);
        x = x(:,end); %ltx only use the last column
    case('fh');     x = huniform(Arg);
    case('pfix');   x = [BdBox([1 3]); BdBox([1 4]); ...
            BdBox([2 3]); BdBox([2 4])];
        
        %ltx for use with PolyMesher and DistMesh
    case('BdBox');  x = BdBox;
    case('MAT');     x = struct('D',IsoElasMtrx(para.E, ...
            para.pos), 'den', para.den); %ltx material
    case('BCond');   x = BoundaryCond(Arg, para, BdBox);
    case('Output');  x = OutputRequests(Arg, BdBox);
    case('EXACT');   x = ExactSln(Arg, para);%ltx exact solution of the problem
    otherwise
        warning('Unexpected keyword in ProbDefPureBendingBeam.')
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
%ltx displacement constrains (prescribed acceleration in a response
%ltx history analysis). One constrain per row: \texttt{[Node Dir]}
%ltx find nodes at left side of beam
lNodes = find(abs(coord(:,1)-BdBox(1))<eps); %ltx nodes at left edge
llNode = find(abs(coord(lNodes,1)-BdBox(1))<eps & ...
    abs(coord(lNodes,2)-BdBox(3))<eps); %ltx lower left corner
n1 = length(lNodes);
ex = ExactSln(coord(lNodes,:), beam); %ltx exact solution at left side
BC_Disp = [ [lNodes; lNodes(llNode)], [ones(n1,1); 2], ...
    [ex(:,1); ex(llNode,2)]]; %ltx displacement boundary condition

%ltx surface traction at right side of beam
%ltx edges of polygon subdomains
[ MeshEdges ] = meshConnectivity( sdConn );
%ltx centres of edges of subdomains
centres = (coord(MeshEdges(:,1),:)+coord(MeshEdges(:,2),:))/2;
%ltx find edges at the right side of the beam
rEdges = MeshEdges(abs(centres(:,1)-BdBox(2))<eps,:);
%ltx surface traction at nodes
ex1 = ExactSln(coord(rEdges(:,1),:), beam);
ex2 = ExactSln(coord(rEdges(:,2),:), beam);
trac = [ex1(:,3)'; zeros(1,length(rEdges)); ex2(:,3)'; ...
    zeros(1,length(rEdges))]; %surface traction on edges
%ltx nodal forces
F = zeros(2*size(coord,1),1); %ltx initialising force vector
%ltx add nodal forces equivalent to surface traction
F = addSurfTraction(coord, rEdges, trac, F);

%ltx output
x = {BC_Disp, F};
end

%ltx {\bf Output Requests}
function [outputs] =OutputRequests(Arg, BdBox)
coord = Arg{1}; %nodal coordinates
%ltx find the node at upper-right corner
inode = find( abs(coord(:,1)-BdBox(2)) <1.d-5 & ...
    abs(coord(:,2)-BdBox(4)) <1.d-5 );
outputs = struct('DispDOF', 2*inode); %ltx request vertical displacement
end

%ltx \textbf{Exact solution} (see Eq.~\eqref{eq:PureBendingBeamExact})
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
p = [   -0.0000     -0.1606;   -0.0000      0.1625;
    0.0000     -0.5000;    0.0000      0.5000;
    0.4080     -0.5000;    0.4093      0.5000;
    0.5137     -0.0004;    0.7988     -0.5000;
    0.7995      0.5000;    1.0000     -0.0002;
    1.2005      0.5000;    1.2012     -0.5000;
    1.4863     -0.0004;    1.5907      0.5000;
    1.5920     -0.5000;    2.0000     -0.5000;
    2.0000      0.5000;    2.0000      0.1625;
    2.0000     -0.1606 ];

%ltx triangular elements
t = [5       1      3 ;    15      16     19 ;
    7       2      1 ;     1       5      7 ;
    5       8      7 ;    18      17     14 ;
    4       2      6 ;     2       7      6 ;
    10       8     12 ;    10       7      8 ;
    12      15     13 ;    13      10     12 ;
    13      15     19 ;    19      18     13 ;
    13      18     14 ;     9       6      7 ;
    7      10      9 ;    11       9     10 ;
    11      13     14 ;    10      13     11 ];
%ltx output
x = {p, t};
end

%ltx {\bf Rectangular mesh} \label{Code_ProbDefPureBendingBeam_rectmesh}
function [x] = RectMesh(para, beam)
%ltx a rectangular mesh

h0 = para.h0; %ltx desired size of square subdomains
L = beam.L; H = beam.H; %ltx dimensions
nx = ceil(L/h0); ny = ceil(H/h0); %ltx number of subdomains
dx = L/nx; dy = H/ny; %ltx actual size of subdomains
[x, y] = meshgrid((0:nx)*dx, (0:ny)*dy-H/2); %ltx grid of nodes
coord = [x(:), y(:)]; %ltx nodal coordinates

%ltx construct scaled boundary subdomains
[ic, ir] = meshgrid(0:nx-1, 0:ny-1); %ltx grid of subdomains
i0 = (ny+1)*(ic(:))+ir(:)+1; %ltx node at low-left corner of a subdomain
sdConn = cell(nx*ny,1);  %ltx initialising connectivity
for ii = 1:nx*ny %ltx loop over subdomains
    i1 = i0(ii);  i2 = i1+ny+1; % lower left and right nodes
    sdConn{ii}=[i1 i2; i2 i2+1; i2+1 i1+1; i1+1 i1];
end
%ltx scaling centre
[x, y] = meshgrid((0.5:nx)*dx, (0.5:ny)*dy-H/2); %ltx grid of points
sdSC = [x(:), y(:)]; %ltx coordinates of scaling centres

%ltx output
x = struct('coord',coord , 'sdConn',{sdConn}, 'sdSC',sdSC );
end