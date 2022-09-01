function [B,U] = bbspline(p,U,u)

%BBSPLINE   B-spline basis functions.
%
%   [B,U] = BBSPLINE(p,[],u) evaluates the p-th degree Bezier spline basis
%   functions at the free parameter values u.
%
%   [B,U] = BBSPLINE(p,N,u), where N>p is a scalar, evaluates the p-th degree
%   B-spline basis functions with N control vertices and a uniform knot vector
%   at the free parameter values u.
%
%   [B,U] = BBSPLINE(p,U,u), where U is a vector of length > 1, evaluates the
%   p-th degree B-spline basis functions with knot vector U at the free
%   parameter values u.
%
%   p     Spline degree.
%   N     Number of control vertices.
%   U     Knot vector.  If necessary, additional knots equal to 0 and 1 are
%         automatically added to the knot vector so that it starts with p+1
%         zeros and ends with p+1 ones, so that splines constructed from the
%         basis functions B start and end at the first and last control vertex.
%         After this operation, the number of knots is equal to n+p+2, where
%         n is the number of spline intervals and N = n+1 is the number of
%         control vertices.
%   u     Free parameter values in the interval [0,1] where the basis functions
%         are evaluated (M * 1) or (1 * M).
%   B     Spline basis functions (M * N).
%   U     Knot vector.
%
%   The 3-th degree B-spline x(u) with 8 control vertices (X,Y) is obtained as
%   follows:
%
%     p = 3;
%     N = 8;
%     u = [0:0.001:1];
%     B = BBSPLINE(p,N,u);
%     X = [0; 1; 2; 5; 5; 6; 7; 8];
%     Y = [0; 2; 1; 1; 4; 6; 3; 0];
%     x = B*X;
%     y = B*Y;
%     figure;
%     plot(x,y,X,Y,'+');
%     axis('equal');
%
%   The (3,3)-th degree B-spline surface x(u,v) with (nu+1)*(nv+1) = 4*4 control
%   vertices (X,Y,Z) is obtained as follows:
%
%     pu = 3;
%     pv = 3;
%     Nu = 4;
%     Nv = 4;
%     u = [0:0.01:1];
%     v = [0:0.02:1];
%     Bu = bbspline(pu,Nu,u);
%     Bv = bbspline(pv,Nv,v);
%     X = [0 0 0 0;
%          3 3 3 3;
%          6 6 6 6;
%          9 9 9 9];
%     Y = [0 3 6 9;
%          0 3 6 9;
%          0 3 6 9;
%          0 3 6 9];
%     Z = [1 1 1 1;
%          5 4 3 2;
%          7 5 5 6;
%          3 5 7 9];
%     x = Bu*X*Bv';
%     y = Bu*Y*Bv';
%     z = Bu*Z*Bv';
%     figure;
%     surf(x,y,z);
%     hold('on');
%     plot3(X,Y,Z,'or');
%     axis('equal');
%
%   See also BSPLINE, BSPLINE2.

% Mattias Schevenels
% April 2017

if numel(U) == 1
  if U<p+1
    error('Number N of control vertices must be larger than the order p or the spline.');
  end
  U = linspace(0,1,U-p+1);
end

if any(U<0) || any(U>1)
  error('Knots must be in the interval [0,1].');
end

u = u(:);
U = U(:);

nmin = length(find((U~=0)&(U~=1)))+p;
while length(U)<nmin+p+2
  U = [0;U;1];
end

n = length(U)-p-2;     % Number of control vertices minus 1
N = n+1;               % Number of control vertices
m = n+p+1;             % Number of knots minus 1
M = length(u);         % Number of evaluation values

B = zeros(M,N);
if p==0
  for i = 0:n
    B((U(i+1)<=u)&(u<U(i+2)),i+1) = 1;
  end
else
  Mr = [bbspline(p-1,U,u) zeros(M,1)];
  for i = 0:n
    D1 = U(i+p+1)-U(i+1);
    if D1 == 0, C1 = 0; else C1 = 1/D1; end
    D2 = U(i+p+2)-U(i+2);
    if D2 == 0, C2 = 0; else C2 = 1/D2; end
    B(:,i+1) = C1*(u-U(i+1)).*Mr(:,i+1) + C2*(U(i+p+2)-u).*Mr(:,i+2);
  end
end
B(u==1,m-p) = 1;





