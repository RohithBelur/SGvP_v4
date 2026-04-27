% Finite element software.
%
% Main program:
%   mtfem
%
% Calculation of stiffness matrix and nodal forces:
%   gk          :  tangential stiffness matrix
%   gfext       :  external nodal forces
%   gfint       :  internal nodal forces
%   gconv       :  compute convergence parameter
%
% Generation of interpolation data:
%   gintp       :  integration point coordinates and weighting factors
%   gn          :  shape functions
%   gdndxi      :  derivatives of shape functions
%
% Data handling:
%   pdat        :  read and process data file
%   gdof        :  allocate global degrees of freedom
%   phist       :  process deformation history array
%   goutp       :  generate output files
