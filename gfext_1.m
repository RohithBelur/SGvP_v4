function fext = gfext_1(elem, nelem, node, eltp, geom, matr, ...
                              data, dof, ndof,  hist, u, du, dt, incriT, DincriT, Dincri, nodf,nnode)
%
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
fext = zeros(ndof, 1);

nnodf = size(nodf,1);
%
% build external temperature nodal forces in element loop
for ielem = 1:nelem
%
%   element type
  eltpe  = eltp(elem(ielem, 1), :);
%
%   global node and dof numbers
  nnodee = feval(eltpe, 'nnodee');
  ndofne = feval(eltpe, 'ndofne');
  inode  = elem(ielem, 4+(1:nnodee));
  dofe   = dof(inode, :)';
  maxndofn = size(dofe, 1);
  idofn  = (1:maxndofn)';
  if all(all(ndofne<=1))
    idof = dofe(ndofne==1)';
  else
    idof   = dofe(idofn(:, ones(1, nnodee)) <= ...
             ndofne(ones(1, maxndofn), :))';
  end
 
%
%   prepare input for element-routine:
  nodee  = node(inode, :);
  geome  = geom(elem(ielem, 2), :);
  matre  = matr(elem(ielem, 3), :);
  datae  = data(elem(ielem, 4), :);
  histe  = hist(ielem, :);
  ue     = u(idof);
  due    = du(idof);
  
%
%   get element internal nodal forces
if strcmp(eltpe(1:2),'pl') && strcmp(eltpe(3:4),'cs')
    ft = zeros(12,1);
else
    ft = feval(eltpe, 'ft', nodee, geome, matre, int64(datae), ...
                         histe, ue, due, dt, incriT(ielem,1), DincriT(ielem,1));
                     
end
%   account for element contribution to external temperature nodal forces
  fext(idof) = fext(idof) + ft;
  
end
  
% apply nodal forces
if nnodf > 0
  fext(dof(nodf(:, 1)+nnode*(nodf(:, 2)-1))) = fext(dof(nodf(:, 1)+nnode*(nodf(:, 2)-1)))+ Dincri * nodf(:, 3);
end
%
%
%
