function [ coord, sdConn, sdSC ] = createSBFEMesh(probdef, ...
                                                  mesher, para)
%Generate polygon S-element mesh
%
%Inputs:
%  probdef - handle of function defining the problem
%  mesher    = 1, DistMesh
%            = 2, PolyMesher
%            = 3, direct input of triangular mesh
%            = 4, direct input of polygon mesh
%            = 'keyword', access user-defined function to input 
%                 or generate S-Element mesh
%  para    - parameters of element size
%            para.h0: Initial edge length (DistMesh only)
%            para.nPolygon: number of polygons (PolyMesher only)
%            nDiv:  number of 2-node line elements per edge
%
%Outputs:
%  coord(i,:)  -  coordinates of node i
%  sdConn{isd,:}(ie,:)  - S-element conncetivity. The nodes of
%                         line element ie in S-element isd.
%  sdSC(isd,:)  - coordinates of scaling centre of S-element isd

switch(mesher)
    
    case (1) %ltx use DistMesh\label{createSBFEMeshCaseDistMesh}
        fd=@(p) probdef('Dist_DistMesh',p); %ltx distance function
        fh=@(p) probdef('fh',p); %ltx scaled edge length
        h0   = para.h0; %ltx initial edge length
        BdBox = probdef('BdBox'); %ltx bounding box
        BdBox = [BdBox(1), BdBox(3); BdBox(2) BdBox(4)];
        pfix = probdef('pfix'); %ltx hard points
        [p, t] = distmesh2d(fd,fh, h0, BdBox, pfix);
        % %         fprintf('%10.4f  %10.4f;\n',p')
        % %         fprintf('%6d  %6d %6d;\n',t')
        %ltx convert a triangular mesh to a polygon mesh
        [ coord,  sdConn, sdSC ] = triToSBFEMesh( p, t );
        
    case (2) %ltx use PolyMesher \label{createSBFEMeshCasePolyMesher}
        %ltx creat polygon mesh
        nPolygon = para.nPolygon; %ltx number of polygons in mesh
        [coord, polygon] = PolyMesher(probdef,nPolygon,100);
        %ltx convert polygons into S-elements
        [sdConn, sdSC] = polygonToSBFEMesh(coord, polygon);
        
    case (3) %ltx use triangular mesh in problem definition file\label{createSBFEMeshCaseTri}
        ArgOut = probdef('TriMesh');
        p = ArgOut{1}; t = ArgOut{2};
        figure
        % %         opt=struct('LineStyle','-', 'LabelEle', 16, ...
        % %             'LabelNode', 14); %ltx plotting options
        PlotTriFEMesh( p, t )
        
        %ltx convert a triangular mesh to a polygon mesh
        [ coord,  sdConn, sdSC ] = triToSBFEMesh( p, t );
        
    case (4) %ltx use polygon mesh in problem definition function\label{createSBFEMeshCasePoly}
        ArgOut = probdef('PolygonMesh');
        coord = ArgOut{1}; polygon = ArgOut{2};
        %ltx convert polygons into S-elements
        [sdConn, sdSC] = polygonToSBFEMesh(coord, polygon);
    otherwise  %ltx obtain mesh from problem definite function\label{createSBFEMeshCaseDef}
        [ glbMesh ] = probdef(mesher, para);
        coord = glbMesh.coord;
        sdConn = glbMesh.sdConn;
        sdSC = glbMesh.sdSC;
        
end

%ltx subdivide one edge into multiple line elements\label{createSBFEMeshSubdiv}
if nargin > 2 && isfield(para,'nDiv') && para.nDiv > 1
    [ meshEdge, sdEdge ] = meshConnectivity( sdConn );
    %ltx cell arrary \texttt{edgeDivXi{i}} contains $\xi$ of edge \texttt{i} 
    %ltx uniform subdivision points in local parametric coordinate $\xi$ (0, 1)
    edgeXi(1:length(meshEdge)) = {( 1:1:para.nDiv-1)'/para.nDiv };
    [coord, sdConn ] = subdivideEdge(edgeXi,...
         coord, meshEdge, sdEdge );
    % %     [ coord, sdConn ] = subdivideAllEdges(...
    % %                            para.nDiv, coord, sdConn );
end

end