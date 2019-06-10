function [ cnt ] = polygonCentroid( xy )
%Centroid of a polygon
%Input
%   xy(i,:) -  coordinates of vertices in counter-clock direction
%Output
%   cnt(:)  -  centroid of polygon

xy = [xy; xy(1,:)]; %ltx appending the 1st vertex after the last one
%ltx array  of all triangles formed by the coordinate origin and an edge
%ltx of the polygon, Eq.~\eqref{eq:triangleArea}
area2 = xy(1:end-1,1).*xy(2:end,2)-xy(2:end,1).*xy(1:end-1,2);
%ltx  array of centroids of triangles, Eq.~\eqref{eq:triangleCentroid}
c = [xy(1:end-1,1)+xy(2:end,1) xy(1:end-1,2)+xy(2:end,2)]/3;
%ltx centroid of polygon, Eq.~\eqref{eq:polygonCentroid}
cnt = sum([area2 area2].*c)/sum(area2);
end