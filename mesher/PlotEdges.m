function  PlotEdges( coord, meshEdge )
%Plot edges
X = [coord(meshEdge(:,1),1)'; coord(meshEdge(:,2),1)'];
Y = [coord(meshEdge(:,1),2)'; coord(meshEdge(:,2),2)'];
plot(X,Y,'b')
meshEdgeCentre = (coord(meshEdge(:,1),:)+coord(meshEdge(:,2),:))/2;
plot(meshEdgeCentre(:,1),meshEdgeCentre(:,2),'x')
nEdge = length(meshEdgeCentre);
text(meshEdgeCentre(:,1),meshEdgeCentre(:,2), ...
      [blanks(nEdge)' int2str((1:nEdge)')],'FontSize',18 );
end

