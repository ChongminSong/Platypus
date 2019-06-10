function [x] = ProbDefLShapedPlate(Demand,Arg)
%Problem definition function of an L-shaped panel
%
%inputs
%   Demand    - option to access data in local functions
%   Arg       - argument to pass to local functions
%output
%   x         - cell array of output. 

%ltx \texttt{L1}: Height of panel; \texttt{L2}: Height of openning;  
%ltx \texttt{E}: Young's modulus; \texttt{pos}: Poisson's ratio; \texttt{den}: mass density    
para = struct('L1',2, 'L2',1, 'E',10E6, 'pos',1/3, 'den',2);
%ltx bounding box $[x_1,x_2,y_1,y_2]$
BdBox = [0 para.L1 0 para.L1]; 
switch(Demand) %ltx access local functions according to \texttt{Demand}
    %ltx direct input triangular mesh
    case('TriMesh'); x = TriMesh; 
        
    %ltx options for mesh generators
    %ltx for use with PolyMesher \label{Code_ProbDefLShapedPlate_PolyMesher} 
    case('Dist');     x = DistFnc(Arg,para); 
    case('BC');       x = {[],[]};
        
    %ltx for use with DistMesh \label{Code_ProbDefLShapedPlate_DistMesh}   
    case('Dist_DistMesh');  x = DistFnc(Arg,para); 
        x = x(:,end); %ltx only use the last column
    case('fh');     x = huniform(Arg);
    case('pfix');  x = [0 para.L1; 0 0; para.L2 para.L2; ...
             para.L2 0; para.L1 para.L2; para.L1 para.L1];     

    %ltx for use with PolyMesher and DistMesh 
    case('BdBox');  x = BdBox;
        
    case('MAT');   x = struct('D',IsoElasMtrx(para.E, ...
            para.pos), 'den', para.den); %ltx material constants
        
    case('BCond');  x = BoundryConditions(Arg,  BdBox);
    case('BCondModal');   x = BoundryConditions(Arg,  BdBox);
    case('BCondTime');   x = BoundryConditions(Arg,  BdBox);
        
    case('Output'); x = OutputRequests(Arg, BdBox);
    case('OutputTime'); x = OutputRequests(Arg, BdBox);
        
    case('EXACT');  x = [];
end
end

%ltx {\bf Compute distance functions}\label{Code_ProbDefLShapedPlate_DistFnc}
function Dist = DistFnc(P,plate)
d1 = dRectangle(P,0,plate.L1,0,plate.L1); %(P,x1,x2,y1,y2)
d2 = dRectangle(P,plate.L2,plate.L1,0,plate.L2);
Dist = dDiff(d1,d2);
end

%ltx {\bf Prescribe boundary conditions}
function [x] = BoundryConditions(Arg,BdBox)
coord = Arg{1}; %ltx nodal coordinates
sdConn = Arg{2}; %ltx element connectivity

%ltx displacement constrains (prescribed acceleration in a response
%ltx  history analysis). One constrain per row: \texttt{[Node Dir]}
eps = 0.001*sqrt((BdBox(2)-BdBox(1))*...
    (BdBox(4)-BdBox(3))/size(coord,1));
rightNodes = find(abs(coord(:,1)-BdBox(2))<eps);
n1 = length(rightNodes);
bottomNodes = find(abs(coord(:,2)-BdBox(3))<eps);
n2 = length(bottomNodes);
% % BC_Disp = [ [rightNodes; bottomNodes; bottomNodes] ...
% %     [ones(n1,1); ones(n2,1); 2*ones(n2,1)] zeros(n1+n2+n2,1)];
BC_Disp = [ [rightNodes; bottomNodes] ...
    [ones(n1,1); 2*ones(n2,1)] zeros(n1+n2,1)];

%ltx surface traction
%ltx edges of polygon S-elements
MeshEdges = meshConnectivity( sdConn );
%ltx centres of edges of S-elements
centres = (coord(MeshEdges(:,1),:)+coord(MeshEdges(:,2),:))/2; 
%ltx edges at the right side of the square
edge = MeshEdges(abs(centres(:,2)-BdBox(4))<eps,:); 
trac = [ 0; 100; 0; 100];  %ltx uniform surface traction
trac = trac(:, ones(1,length(edge))); %ltx surface traction on all edges
%ltx nodal forces
F = zeros(2*size(coord,1),1); %ltx initialising force vector
%ltx add nodal forces equivalent to surface traction 
F = addSurfTraction(coord, edge, trac, F); 
forceHistory = [0 0; 0.001 1; 0.002  0; 200 0]; %ltx see Fig.\ref{fig:L-shaped-panel}b \label{Code_ProbDefLShapedPlate_ft}

%ltx output
x = {BC_Disp, F, forceHistory};
end

%ltx {\bf Output Requests}
function [outputs] =OutputRequests(Arg, BdBox)
coord = Arg{1};
%ltx find the node at upper-right corner
inode = find( abs(coord(:,1)-BdBox(2)) <1.d-5 & ...
        abs(coord(:,2)-BdBox(4)) <1.d-5 ); 
outputs = struct('DispDOF', 2*inode); %ltx request vertical displacement
end

%ltx {\bf Triangular mesh}\label{Code_ProbDefLShapedPlate_TriMesh}
function [x] = TriMesh()
%ltx nodal coordinates
p = [    -0.0000      0.4854;    -0.0000      1.5174;
   -0.0000      1.0308;     0.0000      0.0000;
    0.0000      2.0000;     0.4828      2.0000;
    0.4960      0.7059;     0.4980      0.0000;
    0.5994      1.3949;     0.9679      2.0000;
    1.0000      0.4734;     1.0000      0.0000;
    1.0000      1.0000;     1.2864      1.5042;
    1.5142      2.0000;     1.5267      1.0000;
    2.0000      1.0000;     2.0000      2.0000;
    2.0000      1.5021 ];
%ltx triangles
t = [ 19      14     16 ;    16      14     13 ;
    13       7     11 ;    19      18     15 ;
    15      14     19 ;    13      14      9 ;
     2       3      9 ;     9       7     13 ;
     9       3      7 ;     7       3      1 ;
    19      16     17 ;     8      11      7 ;
     8       1      4 ;     7       1      8 ;
    14      15     10 ;    10       9     14 ;
    11       8     12 ;     9      10      6 ;
     6       5      2 ;     2       9      6 ];
%ltx output
x = {p, t}; 
end