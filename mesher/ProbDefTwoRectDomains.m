function [x] = ProbDefTwoRectDomains(Demand,Arg)
%Problem definition function of a L-shaped panel
%
%inputs
%   Demand    - option to access data in local functions
%   Arg       - argument to pass to local functions
%output
%   x         - cell array of output. 

%ltx \texttt{E}: Young's modulus; \texttt{pos}: Poisson's ratio; \texttt{den}: mass density    
para = struct( 'E',10E6, 'pos',1/3, 'den',2);
%ltx bounding box [x1,x2,y1,y2]
switch(Demand) %ltx access local functions according to \texttt{Demand}
    case('MESH');   x = mesh(Arg);
    case('MAT');   x = struct('D',IsoElasMtrx(para.E, ...
            para.pos), 'den', para.den); %ltx material 
    case('BCond');  x = BoundryConditions(Arg);
    case('BCondModal');  x = BoundryConditions(Arg);
    case('Output'); x = [];
    case('EXACT');  x = [];
end
end

%ltx {\bf Mesh composed of two rectangular domains}
function [x] = mesh (Arg) 

d = cell(2,1); %ltx initilisation of variable storing domains

%ltx mesh of lower rectangular domain
para = Arg{1};
[ d{1}.coord, d{1}.sdConn, d{1}.sdSC ] = ...
    createSBFEMesh(@ProbDefPureBendingBeam, para.mesher, para);
[ d{1}.meshEdge, d{1}.sdEdge] = meshConnectivity( d{1}.sdConn );

%ltx mesh of lower rectangular domain
para = Arg{2};
[ d{2}.coord, d{2}.sdConn, d{2}.sdSC ] = ...
    createSBFEMesh(@ProbDefPureBendingBeam, para.mesher, para);
%ltx scaling domain about a point
cscl = [ 0.4 0 -0.5]; %ltx [scalingFactor, x- and y- coordiantes of rotation centre]
[ d{2}.coord, d{2}.sdSC ] = scaleDomain( d{2}.coord, d{2}.sdSC, cscl );
%ltx shift domain
cshft = [ 0.5 1]; %ltx amount of shifting of x- and y- coordinates
d{2}.coord = bsxfun(@plus, d{2}.coord, cshft);
d{2}.sdSC = bsxfun(@plus, d{2}.sdSC, cshft);% shift scaling centre
[ d{2}.meshEdge, d{2}.sdEdge] = meshConnectivity( d{2}.sdConn );

%ltx Combining the two rectangular domain
eps = 1d-4; %ltx toleration in coordinates
onLine  = @(xy, eps) find(abs(xy(:,2)-0.5) < eps); %ltx contacting line
[ coord, sdConn, sdSC ] = combineMeshOfTwoDomains( d, onLine, eps );

x = {coord, sdConn, sdSC}; %ltx output
end

%ltx {\bf Prescribe boundary conditions}
function [x] = BoundryConditions(Arg)
coord = Arg{1}; %ltx nodal coordinates
sdConn = Arg{2}; %ltx element connectivity

%ltx displacement constrains (prescribed acceleration in a response
%ltx  history analysis). One constrain per row: \texttt{[Node Dir]}
eps = 1d-5;
lNodes = find(abs(coord(:,1)-min(coord(:,1)))<eps);
bNodes = find(abs(coord(:,2)-min(coord(:,2)))<eps);
BC_Disp = [ [bNodes; bNodes], [ones(length(bNodes),1);...
    2*ones(length(bNodes),1)], 0*[bNodes; bNodes] ]; 

%ltx surface traction
%ltx edges of polygon subdomains
[ MeshEdges ] = meshConnectivity( sdConn );
%ltx find edges at the top
ymax = max(coord(:,2));
tEdges = MeshEdges((abs(coord(MeshEdges(:,1),2)-ymax)+ ...
         abs(coord(MeshEdges(:,2),2)-ymax))<eps,:);
%ltx nodal forces
F = zeros(2*size(coord,1),1); %ltx initialising force vector
%ltx add nodal forces equivalent to surface traction
F = addSurfTraction(coord, tEdges, [0; 1; 0; 1], F);

%ltx output
x = {BC_Disp, F};
end