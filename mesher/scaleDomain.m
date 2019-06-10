function [ coord, sdSC ] = scaleDomain( coord, sdSC, cscl )
%Scale a the coordinates of a domain
%
%inputs
%   coord  - nodal coordintes
%   sdSC   - coordinates of scaling centres.
%   cscl   - cscl(1)  : scalingFactor
%            cscl(2:3): x- and y- coordiantes of centre of rotation

%ltx scaling about a point 
coord = cscl(1)*bsxfun(@minus,coord, cscl(2:3));
sdSC = cscl(1)*bsxfun(@minus,sdSC, cscl(2:3));
%ltx recover original position after scaling
coord = bsxfun(@plus,coord, cscl(2:3));
sdSC = bsxfun(@plus,sdSC, cscl(2:3));

end

