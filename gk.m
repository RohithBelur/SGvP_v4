function K = gk(elem, nelem, node, eltp, geom, matr, data, dof, ndof, hist, u, du)
%
% function K = gk(elem, nelem, node, eltp, geom, matr, data, dof, ndof, hist, u, du)
%
% Assembly of stiffness matrix.
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
%   ndof     :  number of degrees of freedom
%   hist     :  deformation history array
%   u        :  current solution
%   du       :  rate of current solution
%
% output:
%   K        :  tangential stiffness matrix
%

%
% initialization
maxnnodee = size(elem, 2) - 4;
maxndofn  = size(dof, 2);
maxndof   = nelem * (maxnnodee * maxndofn)^2;
ik = zeros(maxndof, 1);
jk = zeros(maxndof, 1);
k  = zeros(maxndof, 1);
pk = 0;
%
% build stiffness matrix in element loop
for ielem = 1:nelem
%
%   element type
  eltpe = eltp(elem(ielem, 1), :);
%
%   global node and dof numbers
  nnodee = feval(eltpe, 'nnodee');
  ndofne = feval(eltpe, 'ndofne');
  inode  = elem(ielem, 4+(1:nnodee));
  dofe   = dof(inode, :)';
  ndofe = sum(sum(ndofne));
  idofn  = (1:maxndofn)';
  if all(all(ndofne<=1))
    idof = dofe(ndofne==1)';
  else
    idof   = dofe(idofn(:, ones(1, nnodee)) <= ...
             ndofne(ones(1, maxndofn), :))';
  end
  
  if size(idof,1)==1
    idof=idof(:);
  end
%
%   prepare input for element-routine
  nodee = node(inode, :);
  geome = geom(elem(ielem, 2), :);
  matre = matr(elem(ielem, 3), :);
  datae = data(elem(ielem, 4), :);
  histe = hist(ielem, :);
  ue    = u(idof);
  due   = du(idof);
%
%   get element stiffness matrix
if strcmp(eltpe(1:2),'pl') && strcmp(eltpe(3:4),'cs')
  % use nodal displacements
  Ke = feval(eltpe, 'Ke', nodee, geome, matre, int64(datae), histe, ue, due);
else
  % use rate of nodal displacements
  Ke = feval(eltpe, 'Ke', nodee, geome, matre, int64(datae), histe, ue, due, [0], [0], [0]);
end
    
%     if ielem == 19296
%         K334 = Ke;
%         eval(['save ','K334.mat',' K334' ]);
%         pause
%     elseif ielem == 344
%         K344 = Ke;
%         eval(['save ','K344.mat',' K344' ]);
%     end
    
%
%   store element contribution
  iik     = pk + (1:ndofe^2)';
  iidof   = idof(:, ones(1, ndofe));
  ik(iik) = iidof(:);
  jk(iik) = reshape(iidof', 1, ndofe^2)';
  k(iik)  = Ke(:);
%
%   update pointer
  pk = pk + ndofe^2;
end
%
% assemble stiffness matrix
K = sparse(ik(1:pk), jk(1:pk), k(1:pk), ndof, ndof);
%
%
%
