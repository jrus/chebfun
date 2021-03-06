function varargout = quiver(N, axisLims, varargin)
% QUIVER   Draw phase plot diagrams, based on ODEs specified with CHEBOPS
%
% Calling sequence:
%   H = QUIVER(N, AXIS, 'OPT1', VAL1, ...)
%
% Here, the inputs are:
%   N    : A chebop, whose N.op arguments specifies a second order scalar ODE,
%          or a coupled system of two first order ODEs.
%   AXIS : A 4-vector with elements [XMIN XMAX YMIN YMAX] that specify the
%          rectangular region shown on the phase plot.
%
% It is possible to pass the method various option pairs of the form
% 'OPTIONNAME', OPTIONVALUE. The options supported are:
%   'XPTS'      : An integer, specifying the resolution of the x-axis for the
%                 quiver plot. Default value: 10.
%   'YPTS'      : An integer, specifying the resolution of the y-axis for the
%                 quiver plot. Default value: 10.
%   'NORMALIZE' : A Boolean, which determines whether the arrows on the quiver
%                 plot or normalized all to have the same length. Default value:
%                 false.
%   'SCALE':      By default, quiver automatically scales the arrows to fit
%                 within the grid. By passing a SCALE argument S, the arrows are
%                 fitted within the grid and then stretched by S. Use S = 0 to
%                 plot the arrows without the automatic scaling. See the
%                 documentation for the built-in MATLAB QUIVER for more
%                 information.
%   'LINESPEC'    Specifies options for the quiver plot. The LINESPEC argument
%                 can either be a string supported by the MATLAB PLOT command,
%                 such as 'ro', or a cell array consisting of plot option names
%                 and their values, e.g. {'linewidth', 4}.
%
% The optional output is
%   H   : A quivergroup handle.
%
% Note: The CHEBOP QUIVER command works by reformulating higher order problems
% as coupled first order systems, evaluating the resulting first order system at
% grid that should be interpreted as values of u and u', then calling the
% built-ing MATLAB QUIVER method on the results.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Set default values:
scale = 1;
normalize = false;
u0 = [];
xpts = 20;
ypts = 20;
linespec = {};

% Parse VARARGIN, go through all elements:
while ( ~isempty(varargin) )    
    if ( ~ischar(varargin{1}) && ~isnumeric(varargin{2}) )
        error('followpath:inputArgument','Incorrect options input arguments');
    end
    val = varargin{2};
    switch lower(varargin{1})
        case 'normalize'
            normalize = val;
        case 'xpts'
            xpts = val;
        case 'ypts'
            ypts = val;
        case 'scale'
            scale = val;
        case 'linespec'
            linespec = val;
    end
    
    % Throw away option name and argument and move on:
    varargin(1:2) = [];
end

% Extract the x and y limits:
xl = axisLims(1:2);
yl = axisLims(3:4);

% If ylim is empty, we solve the problem to obtain a range for plotting on:
if ( isempty(xl) )
    u0 = N\0;
    xl = 1.1*minandmax(u0);
end

if ( isempty(yl) )
    if ( isempty(u0) )
        u0 = N\0;
    end
    yl = 1.1*minandmax(diff(u0));
end

% Convert the operator in N to first order.
firstOrderFun = treeVar.toFirstOrder(N.op, 0, N.domain);

% Vectors for constructing a meshgrid:
y1 = linspace(xl(1), xl(end), xpts);
y2 = linspace(yl(1), yl(end), ypts);

% Get a meshgrid for points of interests in the phase plane.
[x,y] = meshgrid(y1, y2);
u = zeros(size(x));
v = zeros(size(x));

% Phase plane portraits really only make sense for autonomous systems, which
% shouldn't depend on t, hence, we simply take t = 0 for evaluating.
t = 0;

% Check if we got passed a system of too high order:
numEquations = length(strfind(func2str(firstOrderFun),';'));

if ( numEquations > 2 )
    error('CHEBFUN:CHEBOP:quiver:tooHighOrder', ...
        ['The ODE passed to chebop/quiver must either be a scalar, second ' ...
        'order ODE or a system of two first order equations.'])
end

% [TODO]: Could probably vectorize this, with reshapes. For now, just loop.
for i = 1:numel(x)
    res = firstOrderFun(t,[x(i); y(i)]);
    u(i) = res(1);
    v(i) = res(2);
end

if ( normalize )
    % Make all arrows equal length:
    nrm = sqrt(u.^2 + v.^2);
    u = u./nrm;
    v = v./nrm;
end

h = quiver(x, y, u, v, scale,linespec{:});

if ( nargout > 0 )
    varargout{1} = h;
else
    
end
