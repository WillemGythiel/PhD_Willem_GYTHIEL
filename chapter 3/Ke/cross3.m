function [out] = cross3(A,B)
out = 0*A;
out(1) = A(2)*B(3)-A(3)*B(2);
out(2) = A(3)*B(1)-A(1)*B(3);
out(3) = A(1)*B(2)-A(2)*B(1);
end