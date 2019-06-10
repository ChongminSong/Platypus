function Poly2triMeshContour(coord, sdSC, SCStrs, Nodestrs, sdConn)
% Contour subdomain Polygonmesh to Trianglemesh stress
%
%Inputs:
%  coord                  - coordinate of nodes
%  sdSC                   - coordinate of scaling centre
%  SCStrs                 - Stress of scaling centre
%  Nodestrs               - Nodal stress
%                           Nodestrs(:,1) sigma x
%                           Nodestrs(:,2) sigma y
%                           Nodestrs(:,3) sigma xy
%  sdConn                 - subdomain connectivity

%%initialisation of triangle mesh connectivity
n=size(Nodestrs,1);
for i=1:size(sdSC)
    sdConn{i}=[sdConn{i} repmat(n+i,[size(sdConn{i},1),1])];
end

Tri=cell2mat(sdConn);
Nodestrs=[Nodestrs;SCStrs];
coord=[coord;sdSC];

figure  %ltx plot stress contour    Nodestrs(:,1) sigma x
patch('Faces',Tri,'Vertices', coord,...
    'FaceVertexCData',Nodestrs(:,1),'FaceColor','interp','EdgeColor','none');
axis equal; axis off;
axis tight; colorbar;
end