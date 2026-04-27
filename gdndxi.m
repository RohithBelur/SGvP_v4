function dNdxi = gdndxi(iplc, xi)
%
% function dNdxi = gdnxi(iplc, xi)
%
% Calculation of derivatives of shape functions with respect to local
% coordinates.
%
% input:
%   iplc     :  interpolation code
%   xi       :  local coordinates of integration point
%
% output:
%   dNdxi    :  array of derivatives
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
    dNdxi = [ -0.50
               0.50 ];
%
%   3 nodes
  elseif (iplc == 'b03')|(iplc == 'l03')
    dNdxi = [    xi1-0.50
              -2*xi1
                 xi1+0.50 ];
%
%   4 nodes
  elseif (iplc == 'b04')|(iplc == 'l04')
    dNdxi = [  -0.50 + xi1 + (-27*xi1^2 + 2*xi1 + 9)/16
                -2*xi1 + (81*xi1^2 + 14*xi1 - 27)/16
                (-81*xi1^2 - 18*xi1 + 27)/16
                0.50 + xi1 + ( 27*xi1^2 + 2*xi1 - 9)/16 ];
%
%   illegal code
  else
    error(sprintf('Illegal interpolation code : %s.', iplc))
  end
%
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
    dNdxi = [  1   0 
               0   1 
              -1  -1 ];
%
%   6 node
  elseif iplc == 't06'
    dNdxi = [  4 * xi1 - 1         0       
               4 * xi2             4 * xi1      
                   0               4 * xi2 - 1    
              -4 * xi2             4 * xi3 - 4 * xi2
              -4 * xi3 + 1        -4 * xi3 + 1
               4 * xi3 - 4 * xi1  -4 * xi1         ];
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
elseif iplc(1)=='q'
   
%
%   local coordinates of integration point
  xi1 = xi(1);
  xi2 = xi(2);
%
%   4 node
  if iplc == 'q04'
    dNdxi = [ -0.25 * (1-xi2)  -0.25 * (1-xi1)
               0.25 * (1-xi2)  -0.25 * (1+xi1)
               0.25 * (1+xi2)   0.25 * (1+xi1)
              -0.25 * (1+xi2)   0.25 * (1-xi1) ];
%
%   8 node
  elseif iplc == 'q08'
    dNdxi = [  0.25*(1-xi2)*(2*xi1+xi2)   0.25*(1-xi1)*(2*xi2+xi1)
                   -(1-xi2)*xi1          -0.50*(1+xi1)*(1-xi1)
               0.25*(1-xi2)*(2*xi1-xi2)   0.25*(1+xi1)*(2*xi2-xi1)
               0.50*(1+xi2)*(1-xi2)           -(1+xi1)*xi2
               0.25*(1+xi2)*(2*xi1+xi2)   0.25*(1+xi1)*(2*xi2+xi1)
                   -(1+xi2)*xi1           0.50*(1+xi1)*(1-xi1)
               0.25*(1+xi2)*(2*xi1-xi2)   0.25*(1-xi1)*(2*xi2-xi1)
              -0.50*(1+xi2)*(1-xi2)           -(1-xi1)*xi2         ];
%
%   9 node
  elseif iplc == 'q09'
    dNdxi = [
      0.25*(1-2*xi1)*   xi2 *(1-xi2)   0.25*   xi1 *(1-xi1)*(1-2*xi2)
                xi1 *   xi2 *(1-xi2)  -0.50*(1+xi1)*(1-xi1)*(1-2*xi2)
     -0.25*(1+2*xi1)*   xi2 *(1-xi2)  -0.25*(1+xi1)*   xi1 *(1-2*xi2)
      0.50*(1+2*xi1)*(1+xi2)*(1-xi2)       -(1+xi1)*   xi1 *     xi2
      0.25*(1+2*xi1)*(1+xi2)*   xi2    0.25*(1+xi1)*   xi1 *(1+2*xi2)
               -xi1 *(1+xi2)*   xi2    0.50*(1+xi1)*(1-xi1)*(1+2*xi2)
     -0.25*(1-2*xi1)*(1+xi2)*   xi2   -0.25*   xi1 *(1-xi1)*(1+2*xi2)
     -0.50*(1-2*xi1)*(1+xi2)*(1-xi2)           xi1 *(1-xi1)*     xi2
     -2.00*     xi1 *(1+xi2)*(1-xi2)  -2.00*(1+xi1)*(1-xi1)*     xi2  ];
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
