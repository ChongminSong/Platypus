function [ ElasMtrx ] = IsoElasMtrx( E, p )
%Elasticity matrix of a isotropic material
% E: Young's modulus; p: Poisson's ratio

ElasMtrx = E/(1-p^2)*[1 p 0;p 1 0;0 0 (1-p)/2];

end

