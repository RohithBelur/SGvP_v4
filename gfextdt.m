function fext = gfext(nodf, nnodf, incri, nnode, dof, ndof)
%
% NOT USED

% function fext = gfext(nodf, nnodf, incri, nnode, dof, ndof)
%
% Computation of external load vector.
%
% input:
%   nodf     :  nodal forces specification
%   nnodf    :  number of nodal forces
%   incri    :  load factor of current increment
%   nnode    :  number of nodes
%   dof      :  global degree of freedom numbering
%   ndof     :  number of degrees of freedom
%
% output:
%   fext     :  external load vector
%

%
% initialize external nodal forces vector
fext = spalloc(ndof, 1, 100);
%
% apply nodal forces
if nnodf > 0
  fext(dof(nodf(:, 1)+nnode*(nodf(:, 2)-1))) = (incri +1)* pdof(:, 3);
end
%
%
%
