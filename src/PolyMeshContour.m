function PolyMeshContour(Nodestrs, polygon, deformed)
% Contour subdomain Polygonmesh stress
%
%Inputs:

%  Nodestrs               - Nodal stress
%                           Nodestrs(:,1) sigma x
%                           Nodestrs(:,2) sigma y
%                           Nodestrs(:,3) sigma xy
%  polygon{i}             - array of vertices of polygon i
%  deformed(i,:)          - coordinate of node i after deformation

 Npolygon=size(polygon);

 polygon = polygon(1:Npolygon)';                
 MaxNVer = max(cellfun(@numel,polygon));     %Max. num. of vertices in mesh
 PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))]; %Pad cells with NaN
 ElemMat = cellfun(PadWNaN,polygon,'UniformOutput',false);
 ElemMat = vertcat(ElemMat{:});              %Create padded polygon matrix

       patch('Faces',ElemMat,'Vertices', deformed,...
       'FaceVertexCData',Nodestrs(:,1),'FaceColor','interp','EdgeColor','none');   
   
       axis equal; axis off;
       axis tight; colorbar;
         
end