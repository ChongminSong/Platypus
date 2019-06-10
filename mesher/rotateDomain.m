function [ coord, sdSC ] = rotateDomain( coord, sdSC, crot )
%Rotate the coordinates of a domain
%
%inputs
%   coord  - nodal coordintes in degree
%   sdSC   - coordinates of scaling centres.
%   crot   - crot(1)  : angle of rotation
%            crot(2:3): x- and y- coordiantes of centre of rotation

a = crot(1);
rot = [cosd(a) sind(a); -sind(a) cosd(a)];

%ltx scaling about a point 
coord = bsxfun(@minus, coord, crot(2:3))*rot;
sdSC  = bsxfun(@minus, sdSC, crot(2:3))*rot;
%ltx recover original position after rotating
coord = bsxfun(@plus,coord, crot(2:3));
sdSC = bsxfun(@plus,sdSC, crot(2:3));

end

