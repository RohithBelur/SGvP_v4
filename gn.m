function N = gn(iplc, xi)
%
% function N = gn(iplc, xi)
%
% Computation of shape function values.
%
% input:
%   iplc     :  interpolation code
%   xi       :  local coordinates of integration point
%
% output:
%   N        :  array of shape functions
%

%
% One-dimensional
% ---------------
%
if (iplc(1) == 'b')|(iplc(1) == 'l')
%
%   local coordinates of integration point
  xi1 = xi(1);
%
%   2 nodes
  if (iplc == 'b02')|(iplc == 'l02')
    N = [  0.50 * (1-xi1)
           0.50 * (1+xi1) ];
%
%   3 nodes
  elseif (iplc == 'b03')|(iplc == 'l03') 
    N = [ -0.50 *    xi1  * (1-xi1)
                  (1+xi1) * (1-xi1)
           0.50 * (1+xi1) *    xi1  ];
%   4 nodes
  elseif (iplc == 'b04')|(iplc == 'l04')
    N = [ 0.5*(1-xi1) - 0.5*(1-xi1^2) + (-9*xi1^3 + xi1^2 + 9*xi1 - 1)/16
          (1-xi1*xi1) + (27*xi1^3 + 7*xi1^2 - 27*xi1 - 7)/16  
          (-27*xi1^3 - 9*xi1^2 + 27*xi1 + 9)/16   
          0.5*(1+xi1) - 0.5*(1-xi1^2) + (9*xi1^3 + xi1^2 - 9*xi1 - 1)/16 ];;       
%
%   illegal code
  else
    error(sprintf('Illegal interpolation code : %s.', iplc))
  end
%
% Two-dimensional, triangular
% ---------------------------
%
elseif iplc(1) == 't'
%
%   local coordinates of integration point
  xi1 = xi(1);
  xi2 = xi(2);
  xi3 = xi(3);
%
%   3 node
  if iplc == 't03' 
    N = [  xi1
           xi2
           xi3 ];
%
%   6 node
  elseif iplc == 't06'
    N = [      xi1 * (2*xi1-1)
           4 * xi1 *    xi2
               xi2 * (2*xi2-1)
           4 * xi2 *    xi3
               xi3 * (2*xi3-1)
           4 * xi1 *    xi3    ];
%
%   illegal code
  else
    error(sprintf('Illegal interpolation code : %s.', iplc))
  end
%
%
% Two-dimensional, quadrilateral
% ------------------------------
%
elseif iplc(1) == 'q'
%
%   local coordinates of integration point
  xi1 = xi(1);
  xi2 = xi(2);
%
%   4 node
  if iplc == 'q04'
    N = [  0.25 * (1-xi1) * (1-xi2)
           0.25 * (1+xi1) * (1-xi2)
           0.25 * (1+xi1) * (1+xi2)
           0.25 * (1-xi1) * (1+xi2) ];
%
%   8 node
  elseif iplc == 'q08'
    N = [ -0.25 * (1-xi1) * (1-xi2) * (1+xi1+xi2)
           0.50 * (1+xi1) * (1-xi1) * (1-xi2)
          -0.25 * (1+xi1) * (1-xi2) * (1-xi1+xi2)
           0.50 * (1+xi1) * (1+xi2) * (1-xi2)
          -0.25 * (1+xi1) * (1+xi2) * (1-xi1-xi2)
           0.50 * (1+xi1) * (1-xi1) * (1+xi2)
          -0.25 * (1-xi1) * (1+xi2) * (1+xi1-xi2)
           0.50 * (1-xi1) * (1+xi2) * (1-xi2)     ];
%
%   9 node
  elseif iplc == 'q09'
    N = [  0.25 *    xi1  * (1-xi1) *    xi2  * (1-xi2)
          -0.50 * (1+xi1) * (1-xi1) *    xi2  * (1-xi2)
          -0.25 * (1+xi1) *    xi1  *    xi2  * (1-xi2)
           0.50 * (1+xi1) *    xi1  * (1+xi2) * (1-xi2)
           0.25 * (1+xi1) *    xi1  * (1+xi2) *    xi2
           0.50 * (1+xi1) * (1-xi1) * (1+xi2) *    xi2
          -0.25 *    xi1  * (1-xi1) * (1+xi2) *    xi2
          -0.50 *    xi1  * (1-xi1) * (1+xi2) * (1-xi2)
                  (1+xi1) * (1-xi1) * (1+xi2) * (1-xi2) ];
%
%   illegal code
  else
    error(sprintf('Illegal interpolation code : %s.', iplc))
  end
%
%
% Unknown Interpolation type
% --------------------------
%
else
  error(sprintf('Illegal interpolation code : %s.', iplc))
end
%
%
%
