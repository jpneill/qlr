\d .lr

///
// normal equations - solution for beta parameters
neq:{(A$/:inv flip[x]$x)$y}

//TODO: add regularised version of normal equations

///
// sum of squared errors
// sse = (y-Xb)^T(y-Xb)
// @param y - vector
// @param X - matrix
// @param b - vector
sse:{[x;y;b]y$y-:x$b}

///
// residual standard error
// @param X - matrix
// @param y - vector
// @param b - vector
// @return - dict with degrees of freedom and residual standard error
rse:{[x;y;b] `dof`rse!(d;sqrt sse[x;y;b]%d:count[y]-count b)}

///
// Standard error of ceofficients
// @param X - matrix
// @param y - vector
// @param b - vector
// @return errors of coefficients b
sec:{[x;y;b]raze {x where y}'[a;{x=/:x}til count a:sqrt(sse[x;y;b]%count[y]-count b)*inv[flip[x]$x]]}

///
// Coeff of determination (aka R-Squared)
// @param x - matrix
// @param y - vector
// @param b - vector
r2:{[x;y;b]a*a:y cor x$b}

///
// Sum of squares
// @param x - vector
sst:{x$x-:avg x}

///
// Adjusted R-Squares
// @param x - matrix
// @param y - vector
// @param b - vector
ar2:{[x;y;b]1-(count[y]-1)*sse[x;y;b]%(count[y]-count b)*sst[y]}

\d .
