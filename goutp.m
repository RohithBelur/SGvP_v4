function goutp(elem, nelem, node, nnode, eltp, geom, matr, data, ...
               dof, ndof, hist, u, du, outp, noutp, outi, nouti, datf, ...
               iincr, nreac, reacti, logfile)
%
% function goutp(elem, nelem, node, nnode, eltp, geom, matr, data, ...
%                dof, ndof, hist, u, outp, noutp, outi, nouti, datf, ...
%                iincr, nreac, reacti)
%
% Writing of analysis results to output files.
%
% input:
%   elem     :  structure connectivity
%   nelem    :  number of elements in structure
%   node     :  nodal coordinates
%   nnode    :  number of nodes
%   eltp     :  element type groups
%   geom     :  geometry groups
%   matr     :  material property groups
%   alpha    :  expansion matrix
%   data     :  element data groups
%   dof      :  degree of freedom numbering
%   ndof     :  number of degrees of freedom
%   hist     :  deformation history array
%   u        :  solution
%   outp     :  output items
%   noutp    :  number of output items
%   outi     :  output increments
%   nouti    :  number of output increments
%   datf     :  datafile name
%   iincr    :  increment number
%
% output:
%   none
%

%
%
% Get output data
% ---------------
%
% initialize output list
if iincr == 0
  outps = 'elem nelem eltp node nnode outp noutp outi nouti';
else
  outps = [];
end
%
% generate output data arrays item by item
for ioutp = 1:noutp
%
%   name of output item
  outpi = deblank(outp(ioutp, :));
%
%   initialize result array and accounting array
  maxnoutpn = 0;
%   if (strcmp(outpi,'volume')==0) && (strcmp(outpi,'vastrn')==0) && (strcmp(outpi,'vastrs')==0) && (strcmp(outpi,'vaplst')==0) && (strcmp(outpi,'volMod')==0) && (strcmp(outpi,'strnMo')==0) && (strcmp(outpi,'strsMo')==0) && (strcmp(outpi,'sigVMo')==0) && (strcmp(outpi,'vasigV')==0) && (strcmp(outpi,'volTwn')==0) && (strcmp(outpi,'strnTw')==0) && (strcmp(outpi,'strsTw')==0)
  if (strcmp(outpi,'stsnT1')==0) && (strcmp(outpi,'stsnT2')==0) && (strcmp(outpi,'stsnT3')==0) && (strcmp(outpi,'stsnT4')==0)

    eval([outpi, ' = gnan(nnode, maxnoutpn);']) 
    
  else
    %scalaire
    eval([outpi, ' = zeros(1, 1);']) 
  end
  nelemn = zeros(nnode, 1);

%  get element data in element loop
  for ielem = 1:nelem
%
% if ((ioutp==14) && (ielem == 3000))
%     pause
% end

%     element type
    eltpe  = eltp(elem(ielem, 1), :);
%
%     global node and dof numbers
    nnodee = feval(eltpe, 'nnodee');
    ndofne = feval(eltpe, 'ndofne');
    inode  = elem(ielem, 4+(1:nnodee));
    dofe   = dof(inode, :)';
    maxndofn = size(dofe, 1);
    idofn  = (1:maxndofn)';
    idof   = dofe(idofn(:, ones(1, nnodee)) <= ...
             ndofne(ones(1, maxndofn), :))';
%
%     prepare input for element routine:
    nodee  = node(inode, :);
    geome  = geom(elem(ielem, 2), :);
    matre  = matr(elem(ielem, 3), :);
    datae  = data(elem(ielem, 4), :);
    histe  = hist(ielem, :);
    ue     = u(idof);
    due    = du(idof);
%
%     call element function
    eval([outpi, 'e = ',eltpe, '(''', outpi, ''', nodee, geome, ', ...
          'matre, int64(datae), histe, ue, due, 0, 0, 0);'])
     if (isempty(eval([eltpe, '(''', outpi, ''', nodee, geome, ', ...
          'matre, int64(datae), histe, ue, due, 0, 0, 0)'])) == 1)
          eval([outpi, 'e = 0;'])
     end
    
%
%     adjust size of output array if necessary
    eval(['noutpn = size(', outpi, 'e, 2);'])
%     if (noutpn > maxnoutpn) && (strcmp(outpi,'volume')==0) && (strcmp(outpi,'vastrn')==0) && (strcmp(outpi,'vastrs')==0) && (strcmp(outpi,'vaplst')==0) && (strcmp(outpi,'volMod')==0) && (strcmp(outpi,'strnMo')==0) && (strcmp(outpi,'strsMo')==0) && (strcmp(outpi,'sigVMo')==0) && (strcmp(outpi,'vasigV')==0) && (strcmp(outpi,'volTwn')==0) && (strcmp(outpi,'strnTw')==0) && (strcmp(outpi,'strsTw')==0)
    if (noutpn > maxnoutpn) && (strcmp(outpi,'stsnT1')==0) && (strcmp(outpi,'stsnT2')==0) && (strcmp(outpi,'stsnT3')==0) && (strcmp(outpi,'stsnT4')==0)


      eval([outpi, ' = [', outpi, ', gnan(nnode, noutpn-maxnoutpn)];'])
      maxnoutpn = noutpn;
    end
%

%     account for element contribution
%     if (noutpn > 0) && (strcmp(outpi,'volume')==0) && (strcmp(outpi,'vastrn')==0) && (strcmp(outpi,'vastrs')==0) && (strcmp(outpi,'vaplst')==0) && (strcmp(outpi,'volMod')==0) && (strcmp(outpi,'strnMo')==0) && (strcmp(outpi,'strsMo')==0) && (strcmp(outpi,'sigVMo')==0) && (strcmp(outpi,'vasigV')==0) && (strcmp(outpi,'volTwn')==0) && (strcmp(outpi,'strnTw')==0) && (strcmp(outpi,'strsTw')==0)
    if (noutpn > 0) && (strcmp(outpi,'stsnT1')==0) && (strcmp(outpi,'stsnT2')==0) && (strcmp(outpi,'stsnT3')==0) && (strcmp(outpi,'stsnT4')==0)


      inoden = inode(find(nelemn(inode) == 0)');
      eval([outpi, '(inoden, 1:noutpn) = ', ...
            'zeros(length(inoden), noutpn);'])
      eval([outpi, '(inode, 1:noutpn) = ', outpi, ...
            '(inode, 1:noutpn) + ', outpi, 'e;'])
      nelemn(inode) = nelemn(inode) + ones(nnodee, 1);
    end
%     if (strcmp(outpi,'volume')==1) || (strcmp(outpi,'vastrn')==1) || (strcmp(outpi,'vastrs')==1) || (strcmp(outpi,'vaplst')==1) || (strcmp(outpi,'volMod')==1) || (strcmp(outpi,'strnMo')==1) || (strcmp(outpi,'strsMo')==1) || (strcmp(outpi,'sigVMo')==1) || (strcmp(outpi,'vasigV')==1) || (strcmp(outpi,'volTwn')==1) || (strcmp(outpi,'strnTw')==1) || (strcmp(outpi,'strsTw')==1)
    if (strcmp(outpi,'stsnT1')==1) || (strcmp(outpi,'stsnT2')==1) || (strcmp(outpi,'stsnT3')==1) || (strcmp(outpi,'stsnT4')==1)

        
      eval([outpi, ' = ', outpi, ...
            ' + ', outpi, 'e;'])
        
        
    end
    
  end
  
%
%   ignore unreferenced nodes
  inelemn = find(nelemn == 0);
  nelemn(inelemn) = ones(size(inelemn, 1), 1);
%
%   perform nodal averaging
%     if (strcmp(outpi,'volume')==0) && (strcmp(outpi,'vastrn')==0) && (strcmp(outpi,'vastrs')==0) && (strcmp(outpi,'vaplst')==0) && (strcmp(outpi,'volMod')==0) && (strcmp(outpi,'strnMo')==0) && (strcmp(outpi,'strsMo')==0) && (strcmp(outpi,'sigVMo')==0) && (strcmp(outpi,'vasigV')==0) && (strcmp(outpi,'volTwn')==0) && (strcmp(outpi,'strnTw')==0) && (strcmp(outpi,'strsTw')==0)
    if (strcmp(outpi,'stsnT1')==0) && (strcmp(outpi,'stsnT2')==0) && (strcmp(outpi,'stsnT3')==0) && (strcmp(outpi,'stsnT4')==0)
        eval([outpi, ' = full(diag(sparse(1.0 ./ nelemn)) * ', outpi, ');'])
    end

%
%   add output item to save list
  outps = [outps, ' ', outpi];
end

%
if nreac
  outps = [outps, ' reacti'];
end

%
%
% Write output file
% -----------------
%
if ~isempty(outps)
%
%   test increment number
  if (iincr < 0) | (iincr >= 100000000)
    error(['Illegal increment number: ', num2str(iincr), '.'])
  end
%
%   construct output filename

  outpf = [datf, ones(1, 4-length(num2str(iincr)))*'0', num2str(iincr)];
  %outpf = ['gs', ones(1, 4-length(num2str(iincr)))*'0', num2str(iincr)];

%
%   write file
  fprintf(1, '  outpf = %9s\n', outpf)
  fprintf(logfile, '  outpf = %9s\n', outpf);
  eval(['save ', outpf, ' ', outps])
end
%
%
%
