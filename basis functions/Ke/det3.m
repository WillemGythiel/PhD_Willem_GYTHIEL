function [out] = det3(A)
out = A(1)*(A(5)*A(9)-A(6)*A(8))+...
      A(4)*(A(3)*A(8)-A(2)*A(9))+...
      A(7)*(A(2)*A(6)-A(3)*A(5));
end  