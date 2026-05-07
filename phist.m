function hist = phist(hist, elem, nelem, eltp, geom, matr, data)
%
% function hist = phist(hist, elem, nelem, eltp, geom, matr, data)
%
% Initializing and processing of deformation history.
%
% input:
%   hist     :  deformation history array (empty in case of
%               initialization)
%   elem     :  structure connectivity
%   nelem    :  number of elements in structure
%   eltp     :  element type groups
%   geom     :  geometry groups
%   matr     :  material property groups
%   data     :  element data groups
%
% output:
%   hist     :  processed deformation history array
%

%
% in case of initialization
if isempty(hist)
  maxnhiste = 1;
  hist = zeros(nelem, maxnhiste);
  for ielem = 1:nelem
%
%     prepare input for element-routine:
    eltpe = eltp(elem(ielem, 1), :);
    geome = geom(elem(ielem, 2), :);
    matre = matr(elem(ielem, 3), :);
    datae = data(elem(ielem, 4), :);
%
%     get element history
    if strcmp(eltpe(1:2),'pl') && strcmp(eltpe(3:4),'cs')
        histe = feval(eltpe, 'histe', [0], geome, matre, int64(datae), [0], [0], [0]);
    else
        histe = feval(eltpe, 'histe', [0], geome, matre, int64(datae), zeros(1,220), [0], [0], [0]);
    end
%
%     adjust size of hist if necessary
    nhiste = size(histe, 2);
    if nhiste > maxnhiste
      hist = [hist, zeros(nelem, nhiste-maxnhiste)];
      maxnhiste = nhiste;
    end
%
%     store element history in global array
    if nhiste > 0
      hist(ielem, 1:nhiste) = histe;
    end
  end
%
% history processing at end of increment
else
  for ielem = 1:nelem
%
%     prepare input for element-routine:
    eltpe = eltp(elem(ielem, 1), :);
    geome = geom(elem(ielem, 2), :);
    matre = matr(elem(ielem, 3), :);
    datae = data(elem(ielem, 4), :);
    histe = hist(ielem, :);
%
%     get element history
    if strcmp(eltpe(1:2),'pl') && strcmp(eltpe(3:4),'cs')
        histe = feval(eltpe, 'histe', [0], geome, matre, int64(datae), histe, [0], [0]);
    else
        histe = feval(eltpe, 'histe', [0], geome, matre, int64(datae), histe, [0], [0], [0]);
    end
%
%     store element history in global array
    nhiste = size(histe,2);
    hist(ielem, 1:nhiste) = histe;
  end
end
%
%
%



