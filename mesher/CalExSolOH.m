function Sol = CalExSolOH(a,r,theta,mat)

if nargin < 4
    part1 = 1-a^2./r.^2.*(3/2*(cos(2*theta))+cos(4*theta));
    part2 = 3/2*a^4./r.^4.*cos(4*theta);
    Sol.StressInCart(1,:) = part1 + part2;

    part1 = -a^2./r.^2.*(1/2*(cos(2*theta)) - cos(4*theta));
    part2 = -3/2*a^4./r.^4.*cos(4*theta);
    Sol.StressInCart(2,:) = part1 + part2;

    part1 = -a^2./r.^2.*(1/2*(sin(2*theta)) + sin(4*theta));
    part2 = +3/2*a^4./r.^4.*sin(4*theta);
    Sol.StressInCart(3,:) = part1 + part2;
else
    mu = mat(1);
    v = mat(2);
    ka = (3-v)/(1+v);
    p1 = r./a*(ka+1).*cos(theta);
    p2 = 2*a./r.*((1+ka)*cos(theta)+cos(3*theta));
    p3 = -2*a^3./r.^3.*cos(3*theta);
    Sol.disp(1,:) = a/8/mu *(p1+p2+p3);

    p1 = r./a*(ka-3).*sin(theta);
    p2 = 2*a./r.*((1-ka)*sin(theta)+sin(3*theta));
    p3 = -2*a^3./r.^3.*sin(3*theta);
    Sol.disp(2,:) = a/8/mu *(p1+p2+p3);
end
