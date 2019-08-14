function [rowi,coli] = near2(X,Y,xi,yi)
%near2 returns indices of values in X and Y that are close to some x,y point. It's similar to 
%find, for nearest neighbors, on a 2D grid. 
%
%% Syntax
% 
%  ind = near2(X,Y,xi,yi)
%  [rowi,coli] = near2(X,Y,xi,yi)
% 
%% Description
% 
% ind = near2(X,Y,xi,yi) returns a linear index corresponding to values
% in X and Y that are closest to the point given by ( xi, yi ). X
% and Y must be 2D grids of equal size and xi and yi must be scalar. 
%
% [rowi,coli] = near2(X,Y,xi,yi) returns the row and column indices 

%% Examples
% % First create a grid of X and Y data:
% 
% [X,Y] = meshgrid(1:5,10:-1:5)
% 
% % Now we can look for X, Y indices corresponding to the point closest
% to (2.1,9.724): 
% 
% [rowi,coli] = near2(X,Y,2.1,9.724)
%
% % Or if you'd prefer linear indices, 
% 
% ind = near2(X,Y,2.1,9.724)
% 
% % If xi,yi describes some point equidistant between multiple grid points
% % in X and Y, all valid minimum-distance indices are returned: 
% 
% [rowi,coli] = near2(X,Y,2.5,9.1)
% 
%% Author Info
% 
% Chad A. Greene of the University of Texas at Austin's Institute for Geophysics
% wrote this in October of 2014.  You can visit Chad over at his internet
% website, www.chadagreene.com.  
% 
% See also find, interp2, ind2sub, sub2ind, and hypot. 

%% Some input checks:

assert(nargin==4,'near2 requires exactly four inputs: gridded X and Y, and scalar values for xi and yi.')
assert(numel(X)==numel(Y),'X and Y must be the same size.') 
assert(size(X,1)==size(Y,1),'X and Y must be the same size.') 
assert(size(X,2)>1,'X and Y must be 2D grids of equal size.')
assert(isscalar(xi)==1,'xi must be a scalar.') 
assert(isscalar(yi)==1,'yi must be a scalar.') 

% extrapolation warnings: 
if xi<min(X(:)) || xi>max(X(:))
    warning('xi lies outside the range of X.')
end
if yi<min(Y(:)) || yi>max(Y(:))
    warning('yi lies outside the range of Y.')
end

%% Calculate distance between (xi,yi) and all points in (X,Y) 

dst = hypot(X-xi,Y-yi);

%% Return index of shortest distance: 

if nargout<2
    % if only one output, return linear index: 
    rowi = find(dst==min(dst(:))); 
else
    [rowi,coli] = find(dst==min(dst(:))); 
end

end

