function [] = PlotResults2NodeEle(sdSln, sdPstP, D)

figure('Color','white')
hold on
axis equal

fct = 1;

sdIntn = cell(length(sdSln),1);
for isd = 1:length(sdSln); % subdomain number
    nNode = length(sdSln{isd}.node);
    xi = (1:-0.02:0).^2; %radial coordinate
    X = zeros(length(xi), nNode); Y = X; Z = X; 
    C = zeros(length(xi), nNode, 5);
    % displacements and strains at the specified raidal coordinate
    for ii= 1:length(xi)
        [nodexy, dsp, strnNode, GPxy, strnEle] = ...
            SubdomainInDispStrain(xi(ii), sdSln{isd}, sdPstP{isd});
        dsp = reshape(dsp,2,[]);
        deformed = nodexy' + fct*dsp;
        X(ii,:) =  deformed(1,:);
        Y(ii,:) =  deformed(2,:);
        C(ii,:,1:2) = dsp';
        C(ii,:,3  ) = sqrt(sum(dsp.^2));
        strs =  D*strnNode; %nodal stresses
        avgs = (strs(1,:)+strs(2,:))/2;
        rs   = sqrt(((strs(1,:)-strs(2,:))/2).^2+strs(3,:).^2);
        C(ii,:,4:8) = [strs; avgs+rs; avgs-rs]';
        sdIntn{isd} = struct('X',X, 'Y',Y, 'C',C);
    end
end
    %invGray = 0.95*(1-gray);
    
for isd = 1:length(sdSln); % subdomain number
    ne = size(sdSln{isd}.conn,1);
    X = sdIntn{isd}.X;
    Y = sdIntn{isd}.Y;
    Z = zeros(size(X));
    C = sdIntn{isd}.C;
    for ie = 1:ne
        nodes = sdSln{isd}.conn(ie,:);
         surface(X(:,nodes),Y(:,nodes),Z(:,nodes),C(:,nodes,7), ...
             'FaceColor','interp',...
             'EdgeColor','none',  ...
             'FaceLighting','phong');
        
%        contourf(X(:,nodes),Y(:,nodes),C(:,nodes,7), [-4:0.5:4], 'LineStyle','none');
%        contourf(X(:,nodes),Y(:,nodes),C(:,nodes,4),20, 'LineStyle','none');
        
    end
    %daspect([1 1 1]);
    %view(-110, 15);
    % text(xy(:,1), xy(:,2), Z(1,1:end-1)',int2str([1:5]'));
    % plot3(X(1,:), Y(1,:), Z(1,:), '-ko');
    % axis off; xlabel('x'); ylabel('x'); zlabel('N');
    %colormap(invGray)
    
    % figure('Color','white')
    % contourf(X,Y,Z, 10, 'LineStyle','none');
    % hold on
    % % text(1.05*(xy(:,1)-0.02), 1.05*xy(:,2),int2str([1:5]'));
    % % plot(X(1,:), Y(1,:), '-ko');
    % axis off; xlabel('x'); ylabel('x');
    % colormap(invGray)
    % colorbar;
    
end

end