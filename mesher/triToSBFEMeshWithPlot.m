function [ coord,  sdConn, sdSC ] = triToSBFEMeshWithPlot( p, t )
%Convert a triangular mesh to an S-element mesh
%
%Inputs  
%  p(i,:)   - coordinates of node i
%  t(i,:)   - nodal numbers of triangle i
%
%Outputs:
%  coord(i,:)   -  coordinates of node i
%  sdConn{isd,:}(ie,:)  - subdomain conncetivity. The nodes of 
%                         line element ie in subdomain isd. 
%  sdSC(isd,:)  - coordinates of scaling centre of subdomain isd

np = length(p); %ltx number of points
nTri = length(t); %ltx number of triangles
%ltx centriods of triangles will be nodes of polygons.
%ltx triangular element numbers will be the nodal number of polygons.
triCnt = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3; %ltx centroids \label{triToSBFEMeshTriCentroid}

%ltx construct data on mesh connectivity\label{triToSBFEMeshConn}
[ meshEdge, ~, edge2sd, node2Edge, node2sd] = ...
                                   meshConnectivity( t );

%plot edges of triangular mesh
figure ('Position', [100 100 340 340])
X = [p(meshEdge(:,1),1)'; p(meshEdge(:,2),1)'];
Y = [p(meshEdge(:,1),2)'; p(meshEdge(:,2),2)'];
plot(X,Y,'k--', 'LineWidth', 1);
hold on
% % plot(p(:,1),p(:,2),'ko', 'LineWidth',1.5)
%plot centroids of triangles
plot(triCnt(:,1),triCnt(:,2),'k*', 'LineWidth',1)
% text(triCnt(:,1),triCnt(:,2),'*','FontSize',16)
% %plot centres of edges
% meshEdgeCentre = (p(meshEdge(:,1),:)+p(meshEdge(:,2),:))/2;
% plot(meshEdgeCentre(:,1),meshEdgeCentre(:,2),'kx', 'LineWidth',1.5)
% % nEdge = length(meshEdgeCentre);
% % text(meshEdgeCentre(:,1),meshEdgeCentre(:,2), ...
% %       [blanks(nEdge)' int2str((1:nEdge)')],'FontSize',14 );
axis equal; axis off;


%ltx number of subdomains connected to an edge
edgeNsd = cellfun(@length,edge2sd);
%ltx list of boundary edges (connected to 1 subdomain only)
bEdge = find(edgeNsd==1);
%ltx midpoints of boundary edges
bEdgeCentre = (p(meshEdge(bEdge,1),:)+p(meshEdge(bEdge,2),:))/2;
%ltx list of points at boundary
bp = unique(meshEdge(bEdge,:));

%ltx include the points in the middle of boundary edges as nodes of subdomains
nbEdge = length(bEdge);
bEdgeNode = nTri + (1:nbEdge)'; %ltx nodal number
bEdgeIdx(bEdge(:)) = (1:nbEdge)'; %ltx index from edge number

hold on
% % % plot edges
% % X = [p(meshEdge(bEdge,1),1)'; p(meshEdge(bEdge,2),1)'];
% % Y = [p(meshEdge(bEdge,1),2)'; p(meshEdge(bEdge,2),2)'];
% % plot(X,Y,'r')
%ltx  plot points on boundary
plot(p(bp,1), p(bp,2),'ko','LineWidth',1.5)
plot(bEdgeCentre(:,1), bEdgeCentre(:,2),'rx', 'LineWidth',1)

%ltx include the points at boundary as nodes of subdomain mesh
nbp = length(bp);
bNode = nTri + nbEdge + (1:nbp)'; %ltx nodal number
bpIdx(bp) = 1:nbp; %ltx index from point number
%ltx nodal coordinates
coord = [triCnt; bEdgeCentre; p(bp,:)];
%ltx show nodal numbers
nNode = length(coord);
 text(coord(:,1),coord(:,2), ...
     [blanks(nNode)' int2str((1:nNode)')] ,'FontSize',9);

%ltx construct polygon subdomains\label{triToSBFEMeshSubdomain}
sdConn = cell(np,1); %ltx initilising connectivity
sdSC   = zeros(np,2); %ltx initilising scaling centre
for ii = 1:np
    if bpIdx(ii) == 0 %ltx interior point \label{triToSBFEMeshInterior}
        node = node2sd{ii}; %ltx subdomains connected to node \texttt{ii}
        %ltx sort nodes in counterclock direction 
        [ node ] = sortNodes(node, coord(node,:), p(ii,:));
        %ltx scaling centre at current point
        sdSC(ii,:) = p(ii,:);
        %ltx line element connectivity in a subdomain
        sdConn{ii}=[node; node(2:end) node(1)]'; 
    else %ltx boundary point, which can become a node or a scaling centre \label{triToSBFEMeshBoundary}
        be =  bEdgeIdx(node2Edge{ii}); %ltx edges connected to node
        nodee = bEdgeNode( be(be~=0) ); %ltx nodes on boundary edges
        %ltx sort the nodes, except for the one at the current point
        node = [ node2sd{ii}; nodee ];
        [ node ] = sortNodes(node, coord(node,:), p(ii,:) );
        %ltx find the 2 boundary nodes in the node list
        idx1 = find(node==nodee(1));
        idx2 = find(node==nodee(2));
        %ltx maintain counterclock direction and rearrange the nodes as: 
        %ltx [boundary node 1, current point (node), boundary node 2, others]
        if abs(idx1-idx2) > 1
            %ltx the 2 boundary nodes are the 1st and last in the list
            node = [node(end) bNode(bpIdx(ii)) node(1:end-1)]; 
        else
            %ltx the 2 boundary nodes are consecutive on the list
            idx = min(idx1,idx2);
            node = [node(idx) bNode(bpIdx(ii)) ...
                    node(idx+1:end) node(1:idx-1) ];
        end
        %ltx internal angle between two boundary edges
        dxy = diff(coord(node(1:3),:)); %ltx $\Delta_x$, $\Delta_y$ of the 1st 2 edge 
        dl = sqrt(sum(dxy.^2,2)); %ltx length
        dxyn = dxy./[dl dl]; %ltx direction cosin
        %ltx angle between 2 boundary edges
        alpha = real(acosd(sum(dxyn(1,:).*dxyn(2,:))));
        beta = 180 - sign(det(dxyn))*alpha; %ltx internal angle
        if beta < 220  %ltx include current point as a node
            %ltx line element connectivity in a subdomain
            sdConn{ii}=[node; node(2:end) node(1)]';
            %ltx select centroid as scaling centre 
            sdSC(ii,:) = polygonCentroid( coord(node,:) ); 
        else %ltx use current point (concave corner) as a scaling centre
            sdSC(ii,:) = p(ii,:); 
            %ltx line element connectivity in a subdomain
            sdConn{ii}=[node(3:end); node(4:end) node(1)]';
        end
    end
end

%plot typical polygons
for ii = [3, 9, 13];
    X = (reshape(coord(sdConn{ii},1),size(sdConn{ii})))';
    Y = (reshape(coord(sdConn{ii},2),size(sdConn{ii})))';
    plot(X, Y, 'k-','Linewidth',2);
    plot(sdSC(ii,1), sdSC(ii,2),'ko');
    plot(sdSC(ii,1), sdSC(ii,2),'k+')
    text(sdSC(ii,1), sdSC(ii,2), [' ' int2str(ii)],'FontSize',12 );
    
end

%ltx remove unconnected nodes \label{triToSBFEMeshRemove}
a = reshape(vertcat(sdConn{:}),[],1); %ltx all nodes
[c, ~, ~] = unique(a); %ltx unique nodes
i(c) = 1:length(c); %ltx new nodal numbers of the unique nodes
coord = coord(c,:); %ltx update the nodal coordinates accordingly
%ltx update line element connectivity in each subdomain
for ii = 1:length(sdConn)
    sdConn{ii} = reshape(i(sdConn{ii}(:)),[],2);
end
end
%ltx {\bf\texttt{function: sortNodes}} \label{triToSBFEMeshSortNodes}
function [ node ] = sortNodes(node, xy, c) 
%ltx sort nodes in counterclock direction around point \texttt{c} 
xy = bsxfun(@minus, xy, c);
ang = atan2(xy(:,2), xy(:,1)); %ltx angular coordinates
[~, ic] = sort(ang); %ltx sort to increasing angular coordinates
node = (node(ic))'; %ltx rearrange nodes
end