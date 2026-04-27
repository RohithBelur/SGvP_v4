function conv = gconv(fres, fext, ndof, D, ndepn)
%
% function conv = gconv(fres, fext, ndof, D, ndepn)
%
% Computation of convergence parameter.
%
% input:
%   fres     :  residual load vector
%   fext     :  external nodal forces
%   ndof     :  number of degrees of freedom
%   D        :  dependency matrix
%   ndepn    :  number of dependencies
%
% output:
%   conv     :  convergence parameter
%

%
% determine norm of residual and external forces

if ndepn == 0
  normfres = norm(fres);
  normfext = norm(fext);
else
  idofd = 1:ndepn;
  idofi = (ndepn+1):ndof;
  normfres = norm(fres(idofi) + D' * fres(idofd));
  normfext = norm(fext(idofi) + D' * fext(idofd));
end
%
% determine convergence parameter
if normfext > 0
  conv = normfres / normfext;
else
  fprintf(1, '\n  Warning: normfext is zero, 1.0 used instead.\n')
  conv = normfres;
end
%
%
%
