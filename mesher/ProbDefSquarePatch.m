function [x] = ProbDefSquarePatch(Demand,Arg)
%Problem definition of patch test on a square

%ltx \texttt{L}: Length; \texttt{D}: Height; 
%ltx \texttt{E}: Young's modulus; \texttt{pos}: Poisson's ratio; \texttt{den}: mass density    
%ltx \texttt{strsCase}: = 1, $\sigma_{x}=1$; = 2, $\sigma_{y}=1$; = 3, $\tau_{xy}=1$
para = struct('L',1, 'D',1, 'E', 1E4, 'pos',0.2, 'den',2, 'sCase',2);
BdBox = [0 para.L 0 para.D]; %ltx bounding box [x1,x2,y1,y2]
switch(Demand)
    case('para');    x = para; %ltx parameters of problem definitions
    %ltx direct input mesh
    case('TriMesh'); x = TriMesh; %ltx triangular mesh
    case('PolygonMesh'); x = PolygonMesh; %ltx polygon mesh

    %ltx mesh generators    
    %ltx for use with PolyMesher
    case('Dist');   x = DistFnc(Arg,BdBox); 
    case('BC');     x = {[],[]};
    %ltx for use with DistMesh   
    case('DistEnd'); x = DistFnc(Arg,BdBox); 
        x = x(:,end); %ltx only use the last column
    case('fh');      x = huniform(Arg);
    case('pfix');    
        x = [BdBox(1) BdBox(3); BdBox(1) BdBox(4); ...
             BdBox(2) BdBox(3); BdBox(2) BdBox(4)];     
    %ltx for use with both PolyMesher and DistDesh 
    case('BdBox');  x = BdBox;
        
    %ltx material property 
    case('MAT');    x = struct('D', IsoElasMtrx(para.E, ...
                      para.pos), 'den', para.den);
                  
    case('BCond');  x = BoundryConditions(Arg, BdBox, para);
    case('Output'); x = OutputRequests(Arg, BdBox);
    case('EXACT');  x = ExactSln(Arg, para);
    otherwise
        warning('Unexpected keyword in ProbDefSquarePatch.')
end
end

%ltx {\bf Compute distance functions}
function Dist = DistFnc(P,BdBox)
Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
end

%ltx {\bf Prescribe boundary conditions}
function [x] = BoundryConditions(Arg,  BdBox, para)
coord = Arg{1}; %ltx nodal coordinates
sdConn = Arg{2}; %ltx element connectivity

%ltx displacement constrains (prescribed acceleration in a response
%ltx  history analysis). One constrain per row: \texttt{[Node Dir]}
eps = 0.001*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))...
       /length(coord));
%ltx node at lower-left corner
llNode = find( abs(coord(:,2)-BdBox(3))<eps & ...
               abs(coord(:,1)-BdBox(1))<eps ); 
%ltx node at lower-right corner
lrNode = find( abs(coord(:,2)-BdBox(3))<eps & ...
               abs(coord(:,1)-BdBox(2))<eps );
BC_Disp = [ llNode 1 0; llNode 2 0; lrNode 2  0];

%ltx surface traction
%ltx identify edges of mesh
[ meshEdge] = meshConnectivity( sdConn );
meshEdgeCentre = (coord(meshEdge(:,1),:)+coord(meshEdge(:,2),:))/2;
%ltx edges (defined by 2 nodes) at the 4 sides of the square
ledge = meshEdge(abs(meshEdgeCentre(:,1)-BdBox(1))<eps,:);
redge = meshEdge(abs(meshEdgeCentre(:,1)-BdBox(2))<eps,:);
bedge = meshEdge(abs(meshEdgeCentre(:,2)-BdBox(3))<eps,:);
tedge = meshEdge(abs(meshEdgeCentre(:,2)-BdBox(4))<eps,:);
%ltx unit vectors corresponding to the edges at the 4 sides
lunt = ones(1,length(ledge)); runt = ones(1,length(redge));
bunt = ones(1,length(bedge)); tunt = ones(1,length(tedge));

switch(para.sCase) %ltx cases of loading
    case(1); %ltx uniform tension along horizontal direction 
        edge = [ledge; redge];
        trac = [ [-lunt runt]; 0*[lunt runt]; ...
                 [-lunt runt]; 0*[lunt runt] ];
    case(2); %ltx uniform tension along vertical direction
        edge = [bedge; tedge];
        trac = [ 0*[bunt tunt]; [-bunt tunt];...
                 0*[bunt tunt]; [-bunt tunt]];
    case(3); %ltx pure shear
        edge = [bedge; redge; tedge; ledge];
        trac = [ [-bunt 0*runt tunt 0*lunt]; ...
                 [0*bunt runt 0*tunt -lunt];...
                 [-bunt 0*runt tunt 0*lunt]; ...
                 [0*bunt runt 0*tunt -lunt] ];
end

%ltx nodal forces
F = zeros(2*size(coord,1),1); %ltx initialising force vector
%ltx add nodal forces equivalent to surface traction 
F = addSurfTraction(coord, edge, trac, F); 

%ltx output
x = {BC_Disp, F};
end

%ltx {\bf Output Request}
function [outputs] =OutputRequests(Arg, BdBox)
coord = Arg{1};
inode = find( abs(coord(:,1)-BdBox(2)) <1.d-5 & ...
              abs(coord(:,2)-BdBox(4)) <1.d-5 );
dof = reshape([2*inode-1 2*inode]',[],1);
outputs = struct('DispDOF', dof);
end

%ltx {\bf Exact solution}
function [x] = ExactSln(xy, para)
x = xy(:,1);  y = xy(:,2);    unt = ones(length(x),1);
E = para.E;   pos = para.pos;  p = 1;
switch(para.sCase)
    case(1); %ltx constant $\sigma_{x}$, see 
        ux = p*x/E;        uy = -p*pos*y/E;
        sx = p*unt;        sy = 0.*unt;     sxy = 0.*unt;
    case(2);  %ltx constant $\sigma_{y}$, see
        ux = -p*pos*x/E;   uy = p*y/E;
        sx = 0.*unt;       sy = p*unt;   sxy = 0.*unt;
    case(3);  %ltx constant $\tau_{xy}$, see
        ux = p*y/( E/(2*(1+pos)) );  uy = 0.*unt;
        sx = 0.*unt;       sy = p*unt;   sxy = 0.*unt;
end
x = [ux, uy, sx, sy, sxy];
end

%ltx {\bf Triangular mesh}
function [x] = TriMesh()
%ltx nodal coordinates
p = [ 0.00  0.00;   0.00  1.00;   0.35  0.65;
      0.48  0.00;   1.00  0.52;   1.00  0.00;
      1.00  1.00];

% % p = [    0.0000      0.0000;    0.0000      0.5028;
% %     0.0000      1.0000;    0.3210      0.4985;
% %     0.5000      0.0000;    0.5000      1.0000;
% %     0.6790      0.4985;    1.0000      0.0000;
% %     1.0000      0.5028;    1.0000      1.0000 ];
%ltx triangles
t = [ 5  4  6;   3  2  1;   1  4  3;
      4  5  3;   2  3  7;   7  3  5];
% % t = [     8       7      5;     9       7      8;
% %      9      10      7;     7      10      6;
% %      1       5      4;     4       2      1;
% %      4       5      7;     7       6      4;
% %      3       2      4;     4       6      3];
%ltx output
x = {p, t}; 
end

%ltx {\bf Polygon mesh}
function [x] = PolygonMesh()
%ltx nodal coordinates
coord = [ 1.00  1.00;   0.57  1.00;   0.00    1.00;
          0.53  0.75;   0.82  0.50;   1.00    0.50;
          1.00  0.00;   0.36  0.63;   0.00    0.77;
          0.53  0.21;   0.33  0.38;   0.55    0.00;
          0.00  0.28;   0.00  0.00];
%ltx polygons
polygon = { [  9   8   4   2   3]; [ 11  10   5   4   8];
            [ 10  12   7   6   5]; [ 14  12  10  11  13];
            [  4   5   6   1   2]; [ 13  11   8   9]    };

%ltx output
x = {coord, polygon}; 

end
