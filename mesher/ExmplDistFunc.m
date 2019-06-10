close all


[X,Y] = meshgrid(0:.05:11, 0:.05:8);
P = [X(:) Y(:)];
b = 1;
 d1 = dCircle(P, 4, 4, 3); %xc,yc,r
 d2 = dRectangle(P, 4, 10, 1, 5);
Dist = d1;
 % d3 = dCircle(P, -1, -1, 0.8); %xc,yc,r
% d4 = dCircle(P, 1, 1, 0.5); %xc,yc,r
% %dLine
% Dist = dUnion( dUnion(d1,d2),d3);
% Dist = dDiff(Dist, d4);
% Dist = dLine(P,6,2,9,8);
% Dist = d2;
% Dist = dUnion(d1,d2);
% Dist = dIntersect(d1,d2);
% Dist = dDiff(d1,d2);
% Dist = dDiff(d2,d1);
Z = reshape(Dist(:,end),size(X,1),[]);
figure('Position', [ 100 100 400 320]);
surfc(X,Y,Z, ...
   'EdgeColor','none')
axis tight; hold on
view(-50,30)
colormap(jet)
colorbar
% surfc(X,Y, Z )
figure('Position', [ 100 100 400 340]);
contour(X,Y, Z , [-2:1:4], 'ShowText','on','LineWidth', 2)
colormap(jet)
%colorbar
axis equal;


Z(Z>1d-1) = NaN;
% figure
% surfc(X,Y,Z,'FaceColor','interp',...
%    'EdgeColor','none')
% axis tight;hold on 
% colorbar
figure('Position', [ 100 100 400 340]);
contour(X,Y, Z , [-3:1:0], 'ShowText','on', 'LineWidth', 2)
colormap(jet)
axis equal; 
