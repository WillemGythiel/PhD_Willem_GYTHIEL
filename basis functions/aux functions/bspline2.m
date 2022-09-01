function [x,dxdX,U,V] = bspline2(X,pu,pv,u,v,U,V)

%BSPLINE2   B-spline surfaces.
%
%   [x,dxdX,U,V] = BSPLINE2(X,pu,pv,u,v,U,V) returns a p-th degree Bezier spline
%   surface with knot vectors (U,V) and control vertices X evaluated at the free
%   parameter values (u,v).
%
%   X     Control vertex coordinates (Nu * Nv).
%   pu    Spline surface degree in the first dimension.
%   pv    Spline surface degree in the second dimension.
%   u     Free parameter values (in the interval [0,1]) in the first dimension
%         used to generate the spline (Mu * 1) or (1 * Mu).
%   v     Free parameter values (in the interval [0,1]) in the second dimension
%         used to generate the spline (Mv * 1) or (1 * Mv).
%   U     Knot vector in the first dimension. See BBSPLINE for more info.
%         Default: uniform knot vector.
%   V     Knot vector in the second dimension. See BBSPLINE for more info.
%         Default: uniform knot vector.
%   x     Spline coordinates (Mu * Mv).
%   dxdX  Derivative of x(:) w.r.t. X(:) ((Mu*Mv) * N).
%   U     Knot vector in the first dimension.
%   V     Knot vector in the second dimension.
%
%   See also BSPLINE, BBSPLINE.

% Mattias Schevenels
% March 2020

Nu = size(X,1);
Nv = size(X,2);

if nargin<6
  U = linspace(0,1,Nu-pu+1);
end
if nargin<7
  V = linspace(0,1,Nv-pv+1);
end

Bu = bbspline(pu,U,u);
Bv = bbspline(pv,V,v);
x = Bu*X*Bv';
dxdX = kron(Bv,Bu);
