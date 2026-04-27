function [pdof, npdof, dof, ndof, D, depn, ndepn, cx1] = gdof(elem, nelem, eltp, nnode, pdof1, npdof1, depn1, ndepn1, iincr, phase_vec, ntransO, cx_p, adof, u_dof, fact)
%
% function [dof, ndof, D] = gdof(elem, nelem, eltp, nnode, pdof, npdof, depn, ndepn)
%
% Allocation of global degrees of freedom.
%
% input:
%   elem     :  structure connectivity
%   nelem    :  number of elements in structure
%   eltp     :  element type groups
%   nnode    :  number of nodes in structure
%   pdof     :  prescribed degrees of freedom
%   npdof    :  number of prescribed degrees of freedom
%   depn     :  dependency specification
%   ndepn    :  number of dependencies
%
% output:
%   dof      :  global position of nodal degrees of freedom
%   ndof     :  number of global degrees of freedom
%   D        :  dependency matrix
%

%
% initialization
ndof = 0;
maxndofn = 1;
dof  = zeros(nnode, maxndofn);
%
%
% Dependent degrees of freedom
% ----------------------------
%
if ndepn1 > 0
%
    ax = find(phase_vec <= iincr);
    ntransN = size(ax,2);
    bx = max(ax);
    cx = find(depn1(:,12) == bx);
    cx1 = [cx_p
           cx];
    cx1 = cx1(find(cx1(:,1)~=0),:);
    
    if ntransN == 1
        depn = depn1(:,1:11);
    elseif ntransN > ntransO        
        depn = depn1(setdiff(1:size(depn1,1),cx1),1:11);
        fprintf(1,'\n  change dpen \n');
    end
    ndepn = size(depn,1);
%
%   enlarge dof if necessary
  maxidofn = max(depn(:, 2));
  if maxidofn > maxndofn
    dof = [dof, zeros(nnode, maxidofn-maxndofn)];
    maxndofn = maxidofn;
  end
%
%   allocate degrees of freedom
  dof(depn(:, 1) + nnode*(depn(:, 2)-1)) = (1:ndepn)';
  ndof = ndof + ndepn;
else
    depn = [];
    cx1 = 0;
    ndepn = 0;
end
%
%
% Prescribed degrees of freedom
% -----------------------------
%
if npdof1 > 0
%
    aa = find(phase_vec <= iincr);
    ntransN = size(aa,2);
    bb = max(aa);
    cc = find(pdof1(:,4) <= bb);
    
    if ntransN > ntransO
        if exist ('adof')
            pdofT = pdof1;
            ch_list = find(pdof1(:,4) == bb);
            ch_dof = adof(pdof1(ch_list,1),3);
            pdofT(ch_list,3) = u_dof(ch_dof)*fact;
            pdof = pdofT(cc,1:3);
            fprintf(1,'\n  change pdof \n');
        else
            pdofT = pdof1(cc,1:3);
            pdof = pdofT(:,:);
            fprintf(1,'\n  change pdof \n');
        end
    end
    npdof = size(pdof, 1);
%
%   enlarge dof if necessary
  maxidofn = max(pdof(:, 2));
  if maxidofn > maxndofn
    dof = [dof, zeros(nnode, maxidofn-maxndofn)];
    maxndofn = maxidofn;
  end

%   allocate degrees of freedom
  dof(pdof(:, 1) + nnode*(pdof(:, 2)-1)) = ndepn + (1:npdof)';
  ndof = ndof + npdof;
end


% Remaining degrees of freedom
% ----------------------------

% generate dof in element loop
for inode=1:nnode

  nnelem=find(sum(elem(:,5:size(elem,2))'==inode)');
  ndofnod=0;
  for jelem = 1:length(nnelem)
    ielem=nnelem(jelem);
%   get relevant element data
    eltpe = eltp(elem(ielem, 1), :);
    ielmnode=find(elem(ielem,5:size(elem,2))==inode);
    ndofne = feval(eltpe, 'ndofne');
    if size(ndofne,1) == 1
      ndofnode=ndofne(ielmnode);
      ndofnod=max(ndofnod,ndofnode);
      mdofnod = ndofnod;
    else
      ndofnode=sum(ndofne(:,ielmnode));
      ndofnod=max(ndofnod,ndofnode);
      mdofnod = max([ndofnod, size(ndofne,1)]);
    end
  end

%   enlarge dof if necessary
  if mdofnod > maxndofn
    dof = [dof, zeros(nnode, mdofnod-maxndofn)];
    maxndofn = mdofnod;
  end
  for idofn = 1:ndofnod
    if size(ndofne,1) == 1
      jdofn = idofn;
    else
      dofi = find(ndofne(:,ielmnode));
      jdofn = dofi(idofn);
    end
    if dof(inode, jdofn) == 0
      ndof = ndof + 1;
      dof(inode, jdofn) = ndof;
    end
  end
end
%
%
% Processing of dependencies
% --------------------------
%
% return empty matrices if no dependencies
if ndepn == 0
  D  = [];
else
%
%   initialization
  maxndofd = floor((size(depn, 2)-2)/3);
  D = spalloc(ndepn, (ndof-ndepn), ndepn*maxndofd);
%
%   construct matrix
  for idepn = 1:ndepn
    ndofd = nnz(depn(idepn, 3*(1:maxndofd)));
    idofd = dof(depn(idepn, 3*(1:ndofd))+nnode*(depn(idepn, 3*(1:ndofd)+1)-1));
    D(idepn, idofd-ndepn) = depn(idepn, 3*(1:ndofd)+2);
  end
end
%

