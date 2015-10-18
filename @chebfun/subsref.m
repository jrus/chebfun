function varargout = subsref(f, index)
%SUBSREF   CHEBFUN subsref.
% ( )
%   F(X) returns the values of the CHEBFUN F evaluated on the array X. If X
%   falls on a breakpoint of F, the corresponding value from F.POINTVALUES is
%   returned. F(X, 'left') or F(X, 'right') will evaluate F at breakpoints
%   using left- or right-hand limits, respectively. See CHEBFUN/FEVAL for
%   further details. F(:) returns F.
%
%   If F is an array-valued column CHEBFUN, and X is an array of arbitrary
%   shape, then F(X) is the result of horizontally concatenating the outputs
%   [F1(X), F2(X), ...] where F1, F2, ... are the columns of F, and each FN(X)
%   has the same shape as X. Likewise, if F is an array-valued row CHEBFUN,
%   then F(X) is the vertically concatenated [F1(X); F2(X); ...].
%
%   If F is an array-valued column CHEBFUN then F(X, COL) returns the values
%   of the columns specified by the vector COL at the points X, again in
%   blocks concatenated horizontally. If F is an array-valued row CHEBFUN then
%   F(ROW, X) returns the values of the rows specified by the vector ROW at
%   the points X, in blocks concatenated vertically.
%
%   Similarly, F(:, COL) or F(ROW, :), for F a column or row CHEBFUN,
%   respectively, returns a new array-valued CHEBFUN containing only the
%   columns or rows specified in COL or ROW.
%
%   F(G), where G is also a CHEBFUN, computes the composition of F and G. See
%   CHEBFUN/COMPOSE for further details.
%
% .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% { }
%   F{S1, S2} restricts F to the domain [S1, S2] < [F.ENDS(1), F.ENDS(end)] and
%   simplifies the result. See RESTRICT  and SIMPLIFY for further details. Note
%   that F{[S1, S2]} is not supported due to the behaviour of the MATLAB
%   subsref() command.
%
% See also FEVAL, COMPOSE, GET, RESTRICT, SIMPLIFY.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Document for array-valued CHEBFUN objects and quasimatrices.

subscripts = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'

        % If the last subscript is 'left' or 'right', save the value
        % for later and then remove it from the list of subscripts.
        extraArguments = {};
        if ( ~isempty(subscripts) )
            lastSub = subscripts{end};
            if ( ischar(lastSub) && ...
                 any(strcmpi(lastSub, {'left', 'right', '-', '+'})))
                extraArguments = {lastSub};
                subscripts(end) = [];
            end
        end

        numCols = numColumns(f);
        isTransposed = f(1).isTransposed;  % True for row chebfuns.

        % Name the array of evaluation points and the array of quasimatrix
        % rows or columns at which to evaluate. If these are not included
        % among the subscripts, default each to ':'.
        switch length(subscripts)
            case 2
                % For row chebfuns, swap the subscripts.
                if ( isTransposed )
                    [columnIndex, x] = subscripts{:};
                else
                    [x, columnIndex] = subscripts{:};
                end

            case 1
                % For both row and column chebfuns, if we get only one
                % subscript, use it to define the evaluation points.
                x = subscripts{1};
                columnIndex = ':';

            case 0
                x = ':';
                columnIndex = ':';

            otherwise
                error('CHEBFUN:CHEBFUN:subsref:dimensions', ...
                    'Index exceeds CHEBFUN dimensions.')
        end

        % Clean up columnIndex and check it for validity. Make sure the index
        % is either the string ':' or a row vector of positive integers, each
        % no greater than the number of columns in f.
        if ( ischar(columnIndex) && isequal(columnIndex, ':') )
            % Do nothing

        elseif ( ~isnumeric(columnIndex) || ...
                 ~isvector(columnIndex) || ...
                 ~isequal(uint32(real(columnIndex)), columnIndex))
            error('CHEBFUN:CHEBFUN:subsref:badsubscript', ...
                 'Column index must be a vector of positive integers.')

        elseif ( max(columnIndex) > numCols )
            error('CHEBFUN:CHEBFUN:subsref:badsubscript', ...
                'Column index exceeds CHEBFUN dimensions.');

        else
            % Reshape column index to be a row vector:
            columnIndex = columnIndex(:).';
            if ( (numel(columnIndex) == numCols) && ...
                  all(columnIndex == 1:numCols) )
                % If columnIndex is all the columns, call it ':' instead.
                columnIndex = ':';
            end

        end

        %% EVALUATE

        if ( isempty(columnIndex) || isempty(x) )
            out = []

        % If x == ':' then return all of the columns requested by columnIndex
        elseif ( ischar(x) && isequal(x, ':') )
            if ( ischar(columnIndex) && isequal(columnIndex, ':') )
                out = f;
            else
                out = extractColumns(f, columnIndex);
            end

        % If x is a chebfun, compose x with f.
        elseif ( isa(x, 'chebfun') )
            out = compose(x, f);

            % TODO: Handle the case where columnIndex != ':', and figure out
            %       the proper behavior for combining row & column chebfuns.

        % If x is a numeric array, evaluate f at x
        elseif ( isnumeric(x) )
            out = feval(f, x, extraArguments{:});

            % If specific columns were requested, extract them.
            % (If columnIndex isn't numeric, it must be ':', so we're done.)
            if ( isnumeric(columnIndex) )

                % Get m, n, p for x of shape [m n p1 p2 p3 ...],
                % with p1*p2*p3*... = p
                [m n p] = size(x);
                outsize = size(x);

                % For row chebfuns we know out is a block array concatenated
                % vertically; for column chebfuns, concatenated horizontally.
                if (isTransposed)
                    out = reshape(out, m, numCols, n*p);
                    outsize(1) = outsize(1) * numel(columnIndex);
                else
                    out = reshape(out, m*n, numCols, p);
                    outsize(2) = outsize(2) * numel(columnIndex);
                end

                % Select the results corresponding to the requested columns,
                % then reshape to the final output size:
                out = out(:, columnIndex, :);
                out = reshape(out, outsize);

            end

        else
            error('CHEBFUN:CHEBFUN:subsref:nonnumeric',...
                'Cannot evaluate chebfun for non-numeric type.')

        end

        % Recurse on SUBSREF():
        if ( numel(index) > 1 )
            out = subsref(out, index(2:end));
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '.'

        % Call GET() for .PROP access.
        out = get(f, subscripts);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRICT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '{}'

        if ( length(subscripts) == 1 )
            if ( isequal(subscripts{1}, ':') )
                % F{:} returns F:
                out = f;
            else
                error('CHEBFUN:CHEBFUN:subsref:badDomain', ...
                    'Invalid domain syntax.')
            end

        elseif ( size(subscripts, 1) == 1 )
            % F{s1,s2,...,sk} returns RESTRICT(F, [s1,s2,...,sk]):
            x = cat(2, subscripts{:});
            out = restrict(f, x);
            out = simplify(out);

        else
            error('CHEBFUN:CHEBFUN:subsref:dimensions', ...
                'Index exceeds chebfun dimensions.')

        end

    otherwise

        error('CHEBFUN:CHEBFUN:subsref:unexpectedType',...
            ['??? Unexpected index.type of ', index(1).type]);
end

% Convert to a cell:
varargout = {out};

end
