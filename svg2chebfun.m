function outcheb = svg2chebfun(pathstring)
%SVG2CHEBFUN  Convert an SVG path to a chebfun.
%   F = SVG2CHEBFUN(P) converts a standard scalable vector graphics path
%   data string to a chebfun, with one piece for each curve segment specified
%   by the svg path.
%
%   SVG is an XML-based vector graphics format implemented by all modern web
%   browsers and available as an export format in most vector graphics
%   software. See the SVG specification for details about SVG path data
%   attributes: http://www.w3.org/TR/SVG/paths.html
%
%   Example:
%       % Plot a triangle one side of which is a cubic curve:
%       path = 'M 0,0   L 1,0   C .5,.1 .1,.5 0,1   L 0,0';
%       f = svg2chebfun(path);
%       plot(f, 'k'), axis equal
%
%       % Rotate the triangle and plot it along with the original:
%       hold on, plot(exp(1i*pi/6) .* f, 'r')
%
%   Note that by default SVG uses a left-handed coordinate system, with
%   the positive y axis pointed downward, whereas MATLAB by default uses
%   a right-handed coordinate system for plotting; it may be necessary
%   to take the complex conjugate of the output of svg2chebfun for its
%   MATLAB plot to match the original SVG drawing.
%
% See also PAFNUTY

% Written by Jacob Rus, October 2015


% Regular expressions for parsing SVG paths via first splitting into command
% groups, and then splitting each group into numeric arguments.
cmdsplit = '\s*([mMzZlLhHvVcCsSqQtTaA])\s*';
numsplit = [ ...
    '\s*,\s*|'                ... % Split at comma with whitespace, or
    '\s+|'                    ... % split at whitespace, or
    '(?<=[0-9])(?=[-+])|'     ... % split before a sign, or
    '(?<=[.eE][0-9]+)(?=[.])' ... % split before a second decimal point.
    ];

% After these lines, cmds is a cell array of commands, e.g. {'M', 'c', 'c'}
% and allparams is a cell array of vectors of command parameters, which
% might represent multiple commands in a row of the same type.
[cmds, allparams] = regexp(pathstring, cmdsplit, 'match', 'split');
cmds = strtrim(cmds);           % Trim any excess whitespace.
allparams = allparams(2:end);   % Ignore part before first command.
allparams = regexp(allparams, numsplit, 'split', 'emptymatch');

% Starting at the point (0, 0), keep a running tally of where the start point
% should be for the next segment. Loop through all of the command blocks,
% and for each one, loop through its constituent segments. For each one,
% add the appropriate command to segmentcmds, and the full set of
% control point locations including the starting point to segmentcoeffs.
% Also keep a running tally of the number of segments, starting at zero.
startpoint = [0 0];
segmentcmds = {};
segmentcoeffs = {};
segments = 0;
for cmdIdx=1:numel(cmds)
    params = str2double(allparams{cmdIdx});
    nparams = numel(params);
    switch cmds{cmdIdx}
        case 'm'  % Move to:
            startpoint = startpoint + params(1:2);
            for k=3:2:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = startpoint + params(k:k+1);
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'M'
            startpoint = params(1:2);
            for k=3:2:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = params(k:k+1);
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'l'  % Line to:
            for k=1:2:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = startpoint + params(k:k+1);
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'L'
            for k=1:2:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = params(k:k+1);
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'h'  % Horizontal line to:
            for k=1:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = startpoint + [params(k) 0];
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'H'
            for k=1:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = [params(k) startpoint(2)];
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'v'  % Vertical line to:
            for k=1:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = startpoint + [0 params(k)];
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'V'
            for k=1:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = [startpoint(1) params(k)];
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'c'  % Cubic curve to:
            for k=1:6:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'C';
                ctrlpt1 = startpoint + params(k:k+1);
                ctrlpt2 = startpoint + params(k+2:k+3);
                endpoint = startpoint + params(k+4:k+5);
                segmentcoeffs{segments} = [startpoint; ctrlpt1; ctrlpt2; endpoint];
                startpoint = endpoint;
            end
        case 'C'
            for k=1:6:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'C';
                ctrlpt1 = params(k:k+1);
                ctrlpt2 = params(k+2:k+3);
                endpoint = params(k+4:k+5);
                segmentcoeffs{segments} = [startpoint; ctrlpt1; ctrlpt2; endpoint];
                startpoint = endpoint;
            end
        case 's'  % Smooth cubic to:
            for k=1:4:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'C';
                if strcmp(segmentcmds{segments-1}, 'C');
                    prevctrlpt2 = segmentcoeffs{segments-1}(3,:);
                    ctrlpt1 = 2 * startpoint - prevctrlpt2;
                else
                    ctrlpt1 = startpoint;
                end
                ctrlpt2 = startpoint + params(k:k+1);
                endpoint = startpoint + params(k+2:k+3);
                segmentcoeffs{segments} = [startpoint; ctrlpt1; ctrlpt2; endpoint];
                startpoint = endpoint;
            end
        case 'S'
            for k=1:4:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'C';
                if strcmp(segmentcmds{segments-1}, 'C');
                    prevctrlpt2 = segmentcoeffs{segments-1}(3,:);
                    ctrlpt1 = 2 * startpoint - prevctrlpt2;
                else
                    ctrlpt1 = startpoint;
                end
                ctrlpt2 = params(k:k+1);
                endpoint = params(k+2:k+3);
                segmentcoeffs{segments} = [startpoint; ctrlpt1; ctrlpt2; endpoint];
                startpoint = endpoint;
            end
        case 'q'  % Quadratic curve to:
            for k=1:4:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'Q';
                ctrlpt = startpoint + params(k:k+1);
                endpoint = startpoint + params(k+2:k+3);
                segmentcoeffs{segments} = [startpoint; ctrlpt; endpoint];
                startpoint = endpoint;
            end
        case 'Q'
            for k=1:4:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'Q';
                ctrlpt = params(k:k+1);
                endpoint = params(k+2:k+3);
                segmentcoeffs{segments} = [startpoint; ctrlpt; endpoint];
                startpoint = endpoint;
            end
        case 't'  % Smooth quadratic to:
            for k=1:2:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'Q';
                if strcmp(segmentcmds{segments-1}, 'Q');
                    prevctrlpt = segmentcoeffs{segments-1}(2,:);
                    ctrlpt = 2 * startpoint - prevctrlpt;
                else
                    ctrlpt = startpoint;
                end
                endpoint = startpoint + params(k:k+1);
                segmentcoeffs{segments} = [startpoint; ctrlpt; endpoint];
                startpoint = endpoint;
            end
        case 'T'
            for k=1:2:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'Q';
                if strcmp(segmentcmds{segments-1}, 'Q');
                    prevctrlpt = segmentcoeffs{segments-1}(2,:);
                    ctrlpt = 2 * startpoint - prevctrlpt;
                else
                    ctrlpt = startpoint;
                end
                endpoint = params(k:k+1);
                segmentcoeffs{segments} = [startpoint; ctrlpt; endpoint];
                startpoint = endpoint;
            end
        case 'a'  % Elliptical arc to:
            % TODO: Implement elliptical arcs, cf.:
            % http://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes
            for k=1:7:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = startpoint + params(k+5:k+6);
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'A'
            for k=1:7:nparams
                segments = segments + 1;
                segmentcmds{segments} = 'L';
                endpoint = params(k+5:k+6);
                segmentcoeffs{segments} = [startpoint; endpoint];
                startpoint = endpoint;
            end
        case 'z'  % Close path:
            % Do nothing in this case.
        case 'Z'
    end
end

% Put out a warning if any of the SVG path segments was an elliptical arc.
% TODO: remove this once elliptical arcs are supported.
if ( any(strcmpi('a', cmds)) )
    warning(['CHEBFUN:SVG2CHEBFUN:SegmentType elliptical arcs ' ...
        'are currently unsupported, treated as line segments.']);
end

% Matrices for converting from Bernstein to Chebyshev basis:
bernstein2cheb3 = 1/32 * [10 6 6 10; -15 -3 3 15; 6 -6 -6 6; -1 3 -3 1];
bernstein2cheb2 = 1/8 * [3 1 3; -4 0 4; 1 -1 1];
bernstein2cheb1 = 1/2 * [1 1; -1 1];

% Multiply all the Bézier curve segment coefficients by the appropriate
% matrices to convert them to Chebyshev basis. Reuse the segmentcoeffs
% cell array for storing the Chebyshev coefficients.
lines = strcmp(segmentcmds, 'L');
if ( sum(lines) )
    cheblines = bernstein2cheb1 * [segmentcoeffs{lines}];
    segmentcoeffs(lines) = mat2cell( ...
        cheblines(:,1:2:end-1) + 1i * cheblines(:,2:2:end), ...
        [2], ones(1, sum(lines)));
end

quads = strcmp(segmentcmds, 'Q');
if ( sum(quads) )
    chebquads = bernstein2cheb2 * [segmentcoeffs{quads}];
    segmentcoeffs(quads) = mat2cell( ...
        chebquads(:,1:2:end-1) + 1i * chebquads(:,2:2:end), ...
        [3], ones(1, sum(quads)));
end

cubics = strcmp(segmentcmds, 'C');
if ( sum(cubics) )
    chebcubics = bernstein2cheb3 * [segmentcoeffs{cubics}];
    segmentcoeffs(cubics) = mat2cell( ...
        chebcubics(:,1:2:end-1) + 1i * chebcubics(:,2:2:end), ...
        [4], ones(1, sum(cubics)));
end

% It is currently impossible to construct a chebfun by passing the chebfun
% constructor a cell array of coefficient vectors, so for now we construct
% a separate single-fun chebfun for each segment, and then pass a cell array
% of those chebfuns back into the chebfun constructor.
%
% TODO: Fix the chebfun constructor to accept a cell array of coefficient
% vectors and then simplify this code. (Should be a big performance boost.)
segmentchebs = {};
for k=1:segments
    segmentchebs{k} = chebfun(segmentcoeffs{k}, [k-1, k], 'coeffs');
end
outcheb = chebfun(segmentchebs, [0:segments]);

end