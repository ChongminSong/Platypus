function [ lineXi ] = findPointsOnLine(line, pnt, res)
%cbgnltx
%{\bf From a given list the points, find those on a line and determine their parametric
% coordinates ($0< \xi < 1$) }
%\par\vspace{1\baselineskip}\par
%{\bf Inputs:}
%\begin{myDescription}
%\item[line(2*i-1:2*i,:)] coordinates of 2 points defining line \texttt{i}
%\item[xyq(i,:)] coordinates of point \texttt{i} on the list of points
%\item[res] resolution (tolerance) of coordinates
%\end{myDescription}
%\vspace{1\baselineskip}\par
%{\bf Outputs:}
%\begin{myDescription}
%\item[\texttt{lineXi\{i\}}] parametric coordinates ($0< \xi < 1$) of
%  points on edge \texttt{i}.
%\end{myDescription}
%cendltx

nline = length(line)/2; %ltx number of lines
lineXi = cell(nline,1); %ltx initialisation
for ii = 1:nline %ltx loop over lines
    xy = line(2*ii-1:2*ii,:);  %ltx coordinates of the two ends of the line
    %ltx shift coordinate orgin to the first point of the line
    xyq = bsxfun(@minus, pnt, xy(1,:));
    xy = xy(2,:)-xy(1,:);
    leng = sqrt(xy(1)^2+xy(2)^2);
    %ltx areas of triangles formed by points in the list and the line\label{matlabCode_findPointsOnLine_area}
    a2 = abs( xyq(:,1)*xy(2) - xyq(:,2)*xy(1) ); 
    %ltx points on line (the heigth of triangle is equal to zero)
    xyq = xyq( a2/leng < res, :); 
    %ltx parametric coordinates of points on line \label{matlabCode_findPointsOnLine_xi}
    if abs(xy(2)) > abs(xy(1)) %ltx use the coordinate with larger change 
        xi = xyq(:,2)/xy(2);  
    else
        xi = xyq(:,1)/xy(1);
    end
    %ltx select the points with their parametric coordinates between 0 and 1
    lineXi{ii} = sort(xi( (1d-5<xi)&(1-xi)>1d-5) );
end

end