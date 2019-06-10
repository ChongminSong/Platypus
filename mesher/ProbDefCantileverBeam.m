function [x] = ProbDefCantileverBeam(Demand,Arg)

%ltx \texttt{L}: Length; \texttt{D}: Height; 
%ltx \texttt{E}: Young's modulus; \texttt{pos}: Poisson's ratio; \texttt{den}: mass density    
para = struct('L',2, 'D',1, 'E', 1E2, 'pos',0.2, 'den',2);
BdBox = [0 para.L -para.D/2 para.D/2]; %ltx bounding box [x1,x2,y1,y2]
h0 = 0.2; %ltx initial mesh size (distmesh2d)
nPolygon = 160; %ltx number of polygons 
switch(Demand)

    case('TriMesh'); x = TriMesh; 

    %ltx for use with PolyMesher
    case('Dist');    
        x = DistFnc(Arg,BdBox); 
    case('BC');      
        x = {[],[]};
    case('nPolygon'); 
        x = nPolygon;
    
    %ltx for use with distmesh2d   
    case('DistEnd'); 
        x = DistFnc(Arg,BdBox); 
        x = x(:,end); %ltx only use the last column
    case('h0');      
        x = h0;
    case('fh');      
        x = huniform(Arg);
    case('pfix');    
        x = [BdBox(1) BdBox(3); BdBox(1) BdBox(4); ...
             BdBox(2) BdBox(3); BdBox(2) BdBox(4)];     
    
    %ltx for use with PolyMesher and distmesh2d 
    case('BdBox');   
        x = BdBox;
        
    %ltx mesher independent 
    case('MAT');     
        x = struct('D',IsoElasMtrx(para.E, para.pos), 'den', para.den);
    case('BCond');    
        x = BoundryConditions(Arg, para, BdBox);
    case('EXACT');   
        x = ExactSln(Arg, para);
end
end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
end

function [x] = BoundryConditions(Arg, beam, BdBox)
coord = Arg.coord;
sdConn = Arg.sdConn;

eps = 0.01*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(coord,1));

%ltx displacement constrains (prescribed acceleration in a response
%ltx  history analysis). One constrain per row: \texttt{[Node Dir]}
leftNodes = find(abs(coord(:,1)-BdBox(1))<eps);
n1 = length(leftNodes);
ex = ExactSln(coord(leftNodes,:), beam);
BC_Disp = [ [leftNodes; leftNodes], ...
    [ones(n1,1); 2*ones(n1,1)], [ex(:,1); ex(:,2)]];

%ltx surface traction
%ltx identify edges of subdomains
[ meshEdge] = meshConnectivity( sdConn );
meshEdgeCentres = (coord(meshEdge(:,1),:)+coord(meshEdge(:,2),:))/2;
%ltx edges at the right side of the square
edge = meshEdge(abs(meshEdgeCentres(:,1)-BdBox(2))<eps,:); %ltx nodes connected to the edges
ex = ExactSln((coord(edge(:,1),:)+coord(edge(:,2),:))/2, beam);
trac = [zeros(1,length(edge));ex(:,4)'; zeros(1,length(edge));ex(:,4)']; 

%ltx output
x = {BC_Disp, edge, trac};
end

function [x] = ExactSln(xy, beam)

x = xy(:,1);  y = xy(:,2);
L = beam.L; D = beam.D;
E = beam.E; pos = beam.pos;
I = D^3/12;
P = 1;
ux = -P*y/(6*E*I).*( (6*L-3*x).*x + (2+pos)*(y.^2-D^2/4) );
uy = P/(6*E*I).*( 3*pos*y.^2.*(L-x) + (4+5*pos)*D^2/4*x + (3*L-x).*x.^2 );
sx = -P*(L-x).*y/I;
sxy = P/(2*I)*(D^2/4-y.^2);
x= [ux, uy, sx, sxy];
end

function [x] = TriMesh()

p = [   -0.0000     -0.1606;
   -0.0000      0.1625;
    0.0000     -0.5000;
    0.0000      0.5000;
    0.4080     -0.5000;
    0.4093      0.5000;
    0.5137     -0.0004;
    0.7988     -0.5000;
    0.7995      0.5000;
    1.0000     -0.0002;
    1.2005      0.5000;
    1.2012     -0.5000;
    1.4863     -0.0004;
    1.5907      0.5000;
    1.5920     -0.5000;
    2.0000     -0.5000;
    2.0000      0.5000;
    2.0000      0.1625;
    2.0000     -0.1606 ];

t = [5       1      3 ;
    15      16     19 ;
     7       2      1 ;
     1       5      7 ;
     5       8      7 ;
    18      17     14 ;
     4       2      6 ;
     2       7      6 ;
    10       8     12 ;
    10       7      8 ;
    12      15     13 ;
    13      10     12 ;
    13      15     19 ;
    19      18     13 ;
    13      18     14 ;
     9       6      7 ;
     7      10      9 ;
    11       9     10 ;
    11      13     14 ;
    10      13     11 ];

x = {p, t}; 
end

