function [fvp] = gfvp(elem, nelem, node, eltp, geom, matr, ...
                              data, dof, ndof,  hist, u, dt)
%
% function [fint, hist] = gfint(elem, nelem, node, eltp, geom, matr, ...
%                               data, dof, ndof, hist, u)
%
% Computation of the viscoplastic force increment Borg 2006 Eqs. 27, 28 RHS
% terms
%
% input:
%   elem     :  structure connectivity
%   nelem    :  number of elements in structure
%   node     :  nodal coordinates
%   eltp     :  element type groups
%   geom     :  geometry groups
%   matr     :  material property groups
%   data     :  element data groups
%   dof      :  global degree of freedom numbering
%   ndof     :  number of construction degrees of freedom
%   hist     :  deformation history array
%   dt       :  time increment
%
% output:
%   fint     :  internal load vector
%   hist     :  deformation history array
%

%
% initialize internal nodal forces

fvp = zeros(ndof, 1);
%
% build internal nodal forces in element loop
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
  


%
%   get element internal nodal forces
if strcmp(eltpe(1:2),'pl') && strcmp(eltpe(3:4),'cs')
  [fvpe] = feval(eltpe, 'fvp', nodee, geome, matre, int64(datae), ...
                         histe, ue, [0]);            % CP
else
  [fvpe] = feval(eltpe, 'fvp', nodee, geome, matre, int64(datae), ...
                         histe, ue, [0], dt, [0], [0]);        % VSGP
end

%     if ielem == 37925
%         pause
%     end

%
%   account for element contribution to internal nodal forces

  fvp(idof) = fvp(idof) + fvpe;
%
if any(isnan(fvpe))
    fprintf(1, ['\n  ielem = %9i  \n'], ielem)
    eval(['save ', 'fvpe.mat fvpe ;'])
    eval(['save ', 'histe.mat histe ;'])
    error('NaN issues in fvp !!!!!\n')
end


end
%
