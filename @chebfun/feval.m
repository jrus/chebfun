function out = feval(F, x, varargin)
%FEVAL   Evaluate a CHEBFUN.
%   FEVAL(F, X) evaluates a CHEBFUN F at the points in X.  If F is a quasimatrix
%   with columns F1, ..., FN, then the result will be [F1(X), ..., FN(X)], the
%   horizontal concatenation of the results of evaluating each column at the
%   points in X.
%
%   FEVAL(F, 'left'), FEVAL(F, 'start'), and FEVAL(F, '-') return the value of F
%   at the left endpoint of its domain.  FEVAL(F, 'right'), FEVAL(F, 'end'), and
%   FEVAL(F, '+') do the same for the right endpoint.
%
%   FEVAL(F, X, 'left') and FEVAL(F, X, '-') evaluate F at the points in X,
%   using left-hand limits to evaluate F at any breakpoints. FEVAL(F, X,
%   'right') and FEVAL(F, X, '+') do the same but using right-hand limits.
%
%   F(X), F('left'), F(X, 'left'), etc, are equivalent syntaxes.
%
%   Example:
%     f = chebfun(@(x) 1./(1 + 25*x.^2));
%     y = feval(f, linspace(-1, 1, 100).');
%
% See also SUBSREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If F or x is empty, there's nothing to do.
if ( isempty(F) )
    out = [];
    return
elseif ( isempty(x) )
    % Return empty matrix with dimensions of the appropriate size.
    out = zeros(size(x));
    return
end

% Deal with FEVAL(FUNCTION_HANDLE, CHEBFUN) type calls.
if ( isa(F, 'function_handle') )
    out = F(x, varargin{:});
    return
end

%% LEFT / RIGHT / END VALUES:
% Support for feval(f, 'left') and feval(f, 'end'), etc.
if ( ischar(x) )
    if ( any(strcmpi(x, {'left', 'start' ,'-'})) )
        out = F.pointValues(1,:);
    elseif ( any(strcmpi(x, {'right', 'end', '+'})) )
        out = F.pointValues(end,:);
    else
        error('CHEBFUN:CHEBFUN:feval:strInput', ...
            'Unknown input argument "%s".', x);
    end
    if ( F(1).isTransposed )
        out = out.';
    end
    return
end

%% EVALUATE:

% For x of shape [m n p1 p2 p3 ...], with p1*p2*p3*... = p, save the
% values m, n, and p for later:
[m, n, p] = size(x);

% This represents the size of our final output. We will adjust it later.
outsize = size(x);

% Reshape x into a single column, for evaluation:
x = x(:);

% Call the number of columns (or rows) in a quasimatrix or array-valued
% chebfun numCols. Evaluate at the column vector of points, yielding
% an output matrix of size (m*n*p, numCols)
numCols = numel(F);
if ( numCols > 1 )
    % If F is a quasimatrix, loop through and evaluate each column:
    out = zeros(m*n*p, numCols);
    for col = 1:numCols
        out(:,col) = columnFeval(F(col), x, varargin{:});
    end
else
    % If F is not a quasimatrix, evaluate all at once:
    out = columnFeval(F, x, varargin{:});
    numCols = size(out, 2);
end

% Reshape the result such that each column will remain a contiguous
% block in our final output, and the third dimension indexes over the
% results from separate rows or columns of a quasimatrix.
if ( F(1).isTransposed )
    out = reshape(out, m, n*p, numCols);
    outsize(1) = outsize(1) * numCols;
else
    out = reshape(out, m*n, p, numCols);
    outsize(2) = outsize(2) * numCols;
end

% Transpose the 2nd and 3rd dimension to properly interleave results from
% separate quasimatrix columns and then reshape output to final size:
out = reshape(permute(out, [1, 3, 2]), outsize);

end

function out = columnFeval(f, x, varargin)
% Evaluate one column of a quasimatrix, or every column of an array-valued
% chebfun. Ignore the distinction between row and column chebfuns.

%% INITIALISE:
dom = f.domain;
funs = f.funs;
numFuns = numel(funs);

%% LEFT AND RIGHT LIMITS:
% Deal with feval(f, x, 'left') and feval(f, x, 'right'):
leftFlag = 0;
rightFlag = 0;
if ( nargin > 2 )
    lr = varargin{1};
    if ( ischar(lr) )
        leftFlag = any(strcmpi(lr, {'left', '-'}));
        rightFlag = any(strcmpi(lr, {'right', '+'}));
        if ( ~(leftFlag || rightFlag) )
            error('CHEBFUN:CHEBFUN:feval:leftRightChar',...
                'Unknown input argument "%s".', lr);
        end
    else
        error('CHEBFUN:CHEBFUN:feval:leftRight', 'Unknown input argument.');
    end
end

%% VALUES FROM FUNS:
if ( numFuns == 1 )
    % Things are simple when there is only a single FUN:
    out = feval(funs{1}, x, varargin{:});
else
    % For multiple FUNs we must determine which FUN corresponds to each x.

    % Preallocate output matrix:
    out = zeros(numel(x), numColumns(f));

    % Replace the first and last domain entries with +/-inf. (Since we want to
    % use FUN{1} if real(x) < dom(1) and FUN{end} if real(x) > dom(end)).
    domInf = [-inf, dom(2:end-1), inf];

    % Loop over each FUN. For all values of x such that real(x) is in
    % [dom(k) dom(k+1)], evaluate FUN{k}(x).
    xReal = real(x);
    for k = 1:numFuns
        I = ( xReal >= domInf(k) ) & ( xReal < domInf(k+1) );
        if ( any(I) )
            out(I,:) = feval(funs{k}, x(I), varargin{:});
        end
    end

end

%% POINTVALUES:
% If the evaluation point corresponds to a breakpoint, we get the value from
% f.pointValues. However, if the 'left' or 'right' flag is given, we first
% modify the entry in point values to take either the left- or right-sided
% limit, respectively.

[xAtBreaks, domIndices] = ismember(x, dom);
if ( any(xAtBreaks) )
    if ( leftFlag )
        % Note that for leftFlag we use 'rval-local', which corresponds to the
        % function value at the right part of the subdomain.
        pointValues = [f.pointValues(1,:); get(f, 'rval-local')];
    elseif ( rightFlag )
        % Similarly rightFlag uses lval-local.
        pointValues = [get(f, 'lval-local'); f.pointValues(end,:)];
    else
        pointValues = f.pointValues;
    end
    % Set the output values for any values of x at the breakpoints.
    out(xAtBreaks,:) = pointValues(domIndices(xAtBreaks),:);
end


end
