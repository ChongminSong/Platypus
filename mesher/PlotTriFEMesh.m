function  PlotTriFEMesh( p, t, opt)
%Plot mesh of triangular elements
%
%Inputs
%    p(i,:)   - coordinates of node i
%    t(i,:)   - nodal numbers of triangle i
%  plot options:
%    opt=struct('LabelEle', 10, 'LabelNode', 8, ...
%               'LineStyle','-','EdgeColor','k','LineWidth', 1)
%  where the options are:
%    opt.LabelEle  : = 0, do not label element;
%                      otherwise, font size of element number
%    opt.LabelNode : = 0, do not label nodes;
%                      otherwise, font size of element number
%    opt.LineStyle, opt.EdgeColor, opt.LineWidth: options of  
%                      'patch' function with the same name

%ltx default options
LineStyle = '-';   EdgeColor = 'k';  LineWidth = 1;
%ltx use specified options if present
if nargin > 2
    if isfield(opt, 'LineStyle') && ~isempty(opt.LineStyle)
        LineStyle = opt.LineStyle;
    end
    if isfield(opt, 'EdgeColor') && ~isempty(opt.EdgeColor)
        EdgeColor = opt.EdgeColor;
    end
    if isfield(opt, 'LineWidth') && ~isempty(opt.LineWidth)
        LineWidth = opt.LineWidth;
    end
end

%ltx plot triangular mesh
patch('Faces',t,'Vertices',p,'FaceColor','w',...
    'LineStyle',LineStyle,'EdgeColor',EdgeColor, ...
    'LineWidth',LineWidth);
axis equal; axis on; grid off;
hold on;

%ltx apply plot options
if nargin > 2
    if isfield(opt, 'LabelNode') && opt.LabelNode %ltx nodes
        if opt.LabelNode <= 2 %ltx font size of nodal number
            fontsize = 14;
        else
            fontsize = opt.LabelNode;
        end
        %ltx plot nodes
        plot(p(:,1),p(:,2),'ko','LineWidth',LineWidth); 
        np = length(p); %ltx number of nodes
        % % text(p(:,1), p(:,2), [blanks(np)' blanks(np)' int2str((1:np)')],'FontSize',11 );
        text(p(:,1), p(:,2), [blanks(np)' int2str((1:np)')], ...
            'FontSize',fontsize ); %ltx label nodes
    end
    if isfield(opt, 'LabelEle') && opt.LabelEle %ltx elements
        if opt.LabelEle <= 2 %ltx font size of element number
            fontsize = 14;
        else
            fontsize = opt.LabelEle;
        end
        %ltx centroids
        triCnt = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3; 
% %        plot(triCnt(:,1),triCnt(:,2),'r*') %ltx plot centroids
        nTri = length(t); %ltx number of triangular elements
% %         text(triCnt(:,1),triCnt(:,2), ...
% %             [blanks(nTri)' int2str((1:nTri)')],...
% %             'FontSize',fontsize, 'Color','r'); %label elements
        text(triCnt(:,1),triCnt(:,2), ...
            [ int2str((1:nTri)')],...
            'FontSize',fontsize, 'Color','r'); %ltx label elements
    end
end