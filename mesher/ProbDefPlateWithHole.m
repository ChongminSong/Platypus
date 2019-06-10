%cbgnltx
%\leftadjust {\bf Problem definition function of a plate with a circular hole, 
% Example~\ref{Example-PlatWithHole}}  
%cendltx
function [x] = ProbDefPlateWithHole(Demand,Arg)
%inputs
%   Demand    - option to access data in local functions
%   Arg       - argument to pass to local functions
%output
%   x         - cell array of output. 

%ltx \texttt{L}: Length; \texttt{D}: Height; \texttt{R}: Radius of hole; 
%ltx \texttt{E}: Young's modulus; \texttt{pos}: Poisson's ratio; \texttt{den}: mass density    
para = struct('L',2, 'D',2, 'R', 0.5, ...
    'E', 1E2, 'pos',0.25, 'den', 2);
BdBox = [-para.L/2 para.L/2 -para.D/2 para.D/2]; % x1,x2,y1,y2
switch(Demand)

    %ltx direct input triangular mesh
    case('TriMesh'); x = []; 

    %ltx mesh generators
    %ltx for use with PolyMesher
    case('Dist');    
        x = DistFnc(Arg,para,BdBox); 
    case('BC');      
        x = {[],[]};
    case('nPolygon'); 
        x = nPolygon;
    
    %ltx for use with distmesh2d   
    case('Dist_DistMesh'); 
        x = DistFnc(Arg,para,BdBox); 
        x = x(:,end); %ltx only use the last column
    case('fh');      
 % % x = huniform(Arg);
        x = 1+2*dcircle(Arg,0,0,para.R);
    case('pfix');    
        x = [BdBox(1) BdBox(3); BdBox(1) BdBox(4); ...
             BdBox(2) BdBox(3); BdBox(2) BdBox(4)];     
    
    %ltx for use with PolyMesher and distmesh2d 
    case('BdBox');   
        x = BdBox;
        
    %ltx mesher independent 
    case('MAT');     
        x = struct('D',IsoElasMtrx(para.E, para.pos), ...
            'den', para.den);  %ltx material
    case('BCond');    
        x = BoundryConditions(Arg,  BdBox);
    case('Output'); x = OutputRequests(Arg, BdBox);

    case('EXACT');   
        x =  {};
end
end
%ltx \textbf{Signed distance function}\label{Code_ProbDefPlatWithHole_DistFnc}
function Dist = DistFnc(P,plate,BdBox)
d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
d2 = dCircle(P, 0, 0, plate.R); %xc,yc,r
Dist = dDiff(d1,d2);
end
%ltx \textbf{Boundary conditions}  \label{Code_ProbDefPlatWithHole_BC}
function [x] = BoundryConditions(Arg,BdBox)
coord = Arg{1}; %ltx nodal coordinates
sdConn = Arg{2}; %ltx element connectivity

eps = 0.01*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(coord,1));

%ltx displacement constrains (prescribed acceleration in a response
%ltx  history analysis). One constrain per row: \texttt{[Node Dir]}
leftNodes = find(abs(coord(:,1)-BdBox(1))<eps);
leftlowerNode = find(abs(coord(:,1)-BdBox(1))<eps & ...
    abs(coord(:,2)-BdBox(3)) <eps);
n1 = length(leftNodes);
BC_Disp = [ [leftNodes; leftlowerNode] ...
    [ones(n1,1); 2] zeros(n1+1,1)];

%ltx surface traction
%ltx identify edges of subdomains
[ meshEdges] = meshConnectivity( sdConn );
meshEdgeCentres = (coord(meshEdges(:,1),:)+coord(meshEdges(:,2),:))/2;
%ltx edges at the right side of the square
edge = meshEdges(abs(meshEdgeCentres(:,1)-BdBox(2))<eps,:); %ltx 1st nodes of edges
trac = [ 1; 0; 1; 0];  %ltx uniform surface traction
trac = trac(:, ones(1,length(edge))); %ltx surface traction on all edges

%ltx nodal forces
F = zeros(2*size(coord,1),1); %ltx initialising force vector
%ltx add nodal forces equivalent to surface traction 
F = addSurfTraction(coord, edge, trac, F); 

%ltx output
x = {BC_Disp, F};
end

%ltx {\bf Output Requests}
function [outputs] = OutputRequests(Arg, BdBox)
coord = Arg.coord;
inode = find( abs(coord(:,1)-BdBox(2)) <1.d-5 & ...
        abs(coord(:,2)-BdBox(4)) <1.d-5 ); %ltx node at upper-right corner
outputs = struct('DispDOF', [2*inode-1; 2*inode]); %ltx DOFs in one column
end