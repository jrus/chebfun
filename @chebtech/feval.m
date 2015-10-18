function y = feval(f, x)
%FEVAL   Evaluate a CHEBTECH.
%   Y = FEVAL(F, X) Evaluation of the CHEBTECH F at points X via Clenshaw's
%   algorithm (see CLENSHAW.m).
%
%   X must be a vector, and the output of FEVAL will be a column vector of
%   the same length as X. To evaluate a CHEBTECH at an array of points with
%   arbitrary shape, call FEVAL on a high-level CHEBFUN object.
%
%   If size(F, 2) > 1 then FEVAL returns values in the form [F_1(X), F_2(X),
%   ...], where length(F_k(X)) = length(X).
%
%   Example:
%     f = chebtech2(@(x) 1./(1 + 25*x.^2));
%     x = linspace(-1, 1, 1000);
%     fx = feval(f, x);
%     plot(x,fx,'r-')
%
% See also BARY, CLENSHAW.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


if ( ~isvector(x) )
    error('CHEBFUN:CHEBTECH:feval:evalArrayAtNDArray', ...
        ['Evaluation of a CHEBTECH requires a vector of inputs; ' ...
         'to evaluate an arbitrary shape, call feval on a high-level ' ...
         'CHEBFUN object.']);
end

if ( isempty(f) || isempty(x) )
    % Empty array with the proper dimensions:
    y = zeros(length(x), size(f, 2);
else
     y = f.clenshaw(x(:), f.coeffs);
end

end
