function [fint, hist] = gfint(elem, nelem, node, eltp, geom, matr, ...
                              data, dof, ndof,  hist, u, du, dt, incriT, DincriT,iincr)
%
% function [fint, hist] = gfint(elem, nelem, node, eltp, geom, matr, ...
%                               data, dof, ndof, hist, u)
%
% Assembly of internal nodal forces.
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
%   u        :  nodal displacement solution
%   du       :  rate of nodal displacement solution
%   dt       :  time increment
%
% output:
%   fint     :  internal load vector
%   hist     :  deformation history array
%

%
% initialize internal nodal forces
fint = zeros(ndof, 1);
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
  due    = du(idof);
%

% if (ielem == 716) && (iincr == )
%     pause
% end


%   get element internal nodal forces
if strcmp(eltpe(1:2),'pl') && strcmp(eltpe(3:4),'cs')
% use nodal displacements
  [finte, histe] = feval(eltpe, 'finte', nodee, geome, matre, int64(datae), ...
                         histe, ue, due);
else
% use rate of nodal displacements
  [finte, histe] = feval(eltpe, 'finte', nodee, geome, matre, int64(datae), ...
                         histe, ue, due, dt, incriT(ielem,1), DincriT(ielem,1));
end

% if (iincr == 5001) && (ielem == 248)
% %     fprintf(1, '\n  ielem = %9i\n', ielem)
% %     eltpe = 'plbb2t15';
%     G = 0.2/(2*(1+0.3));
%     matre = [0.2   0.3   5   0.04   1e-3     0.439e-3  0.25e-3   G/20  G/20   1.0e-3   1.0e-3   0         0        0       0.04    5];
%     [finte, histe] = feval(eltpe, 'finte', nodee, geome, matre, int64(datae), ...
%                          histe, ue, due, dt, incriT(ielem,1), DincriT(ielem,1));
%     
% %     pause
% end

%
%   account for element contribution to internal nodal forces
  fint(idof) = fint(idof) + finte;
%
%   store new history information
  if size(histe,2)<size(hist,2)
      histe = [histe, zeros(1,size(hist,2)-size(histe,2))];
  end
  hist(ielem, :) = histe;
  
if any(isnan(histe))
    fprintf(1, ['\n  ielem = %9i  \n'], ielem)
    eval(['save ', 'histe.mat histe ;'])
    error('NaN issues again !!!!!\n')
end

end