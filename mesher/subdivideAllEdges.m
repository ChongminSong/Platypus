function [coord,  sdConnNew ] = subdivideAllEdges( ...
                                  nDiv, coord, sdConn )
%Subdivide all edges of subdomains into multiple elements

[ meshEdge, sdEdge ] = meshConnectivity( sdConn );

nEdge0 = length(meshEdge); %ltx number of edges in original mesh
edgeDiv = nDiv*ones(nEdge0,1); %ltx number of division of each edge

newEdges = cell(nEdge0,1); %ltx to store edges after division
ib = length(coord) + 1; %ltx index for  new nodes
nExtraPnt = sum(edgeDiv-1); %ltx number of extra points after subdivision
coord = [coord; zeros(nExtraPnt, 2)]; %ltx including new nodes
for ii = 1:nEdge0 %ltx dividing edges in orginial mesh
    if edgeDiv(ii) > 1
        ie = ib + edgeDiv(ii) - 2; %ltx index for the last new node
        dxi = 1/edgeDiv(ii);     %ltx uniform subdivision
        xi = (dxi : dxi : 1-dxi/2)'; %ltx interpolate between 0 and 1
        xyb = coord(meshEdge(ii,1),:); %ltx coordinates of start point
        xye = coord(meshEdge(ii,2),:); %ltx coordinates of end point
        coord(ib:ie,:) = bsxfun(@times,xyb,1-xi) + ...
            bsxfun(@times,xye,xi); %ltx new nodes by interpolation
        newEdges{ii} = [meshEdge(ii,1), ib:ie; ...
            ib:ie meshEdge(ii,2)]'; %ltx new edges on an old one
        ib = ie + 1;
    else
        newEdges{ii} = meshEdge(ii,:); %ltx no division
    end
end

nsd = length(sdEdge); %ltx number of subdomains
sdConnNew = cell(nsd,1);  %ltx initialising connectivity
for isd = 1:nsd
    sdNewEdges = newEdges(abs(sdEdge{isd})); %ltx update connectivity
    %ltx find edges that follow clockwise direction
    idx = find(sdEdge{isd}(:) < 0);
    for ii = idx'
        %ltx reverse the element orientation
        sdNewEdges{ii} = [sdNewEdges{ii}(end:-1:1,2), ...
            sdNewEdges{ii}(end:-1:1,1)];
    end
    %ltx new element connectivity of a subdomain
    sdConnNew{isd} = vertcat(sdNewEdges{:}); 
end

end