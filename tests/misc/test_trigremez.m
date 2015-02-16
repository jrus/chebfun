% Test for trigBary.m.
function pass = test_trigremez(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Set an error tolerance.
tol = 1.0e-12;

%% check for emptiness:
p = trigremez(chebfun);
pass(1) = isempty(p);


%% check for a simple trigfun 
fh = @(x) cos(4*pi*x) + 1;
f = chebfun(fh, [0, 2], 'trig');
n = 1;
[p, err, status] = trigremez(f, n);
error = f - p;
pass(1) = isPeriodicTech(p) & (length(p) == 2*n+1);
pass(2) = norm(abs(error(status.xk)) - err, inf) < 100*tol;
pass(3) = norm(p-1,inf) < 100*tol;
n = 4;

%% Pass a non-trig chebfun:
f = chebfun(f);
[p, err, status] = trigremez(f, n);
error = f - p;
pass(4) = isPeriodicTech(p) & (length(p) == 2*n+1);
pass(5) = norm(abs(error(status.xk)) - err, inf) < 100*tol;
pass(6) = norm(p-f,inf) < 100*tol;

%% check for a funciton with a kink:
a = 1;
b = 5;
s = a + .7*(b-a);
fh = @(x) (s-x)./(s-a).*(x<s) + (x-s)./(b-s).*(x>=s);
f = chebfun(@(x) fh(x), [a, b], 'splitting', 'on');
n = 2;
[p, err, status] = trigremez(f, n);
error = p-f;
xk = status.xk;
pass(7) = norm(abs(error(xk))-err, inf) < 100*tol;
end