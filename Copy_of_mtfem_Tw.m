function mtfem_Tw(datf)
%

% reaction forces computed from dependencies for one node only


% function mtfem(datf)
%
% MTFEM
%
% Nonlinear structural finite eleMent analysis.
% Eindhoven University of Technology, Mechanical Engineering
%
% input:
%   datf     :  model data file
%
%profile on
%
% Processing of model data
% ------------------------
%
% show name of datafile
fprintf(1, '  datf  = %9s\n', datf)
% create a logfile
logname  = strcat(datf,'.log');
logfile = fopen(logname, 'w');
%
% load model data
[node, nnode, elem, nelem, eltp, geom, matr, data, ...
 pdof1, npdof1, nodf, nnodf, depn1, ndepn1, incr, nincr, ...
 tolr, maxniter, reac, nreac, outp, noutp, outi, ...
 nouti, restart, time, phase_vec, incrTt, nodes_t] = pdat(datf);
%
%
iincr = 0;
ntransO = 0;
cx_p = 0;
%
% allocate nodal degrees of freedom
[pdof, npdof, dof, ndof, D] = gdof(elem, nelem, eltp, nnode, pdof1, npdof1, depn1, ndepn1, ...
                                                             iincr, phase_vec, ntransO, cx_p);

% show some info on the model
fprintf(1, ['\n  nelem = %9i  nnode = %9i  ndof  = %9i', ...
            '  npdof = %9i\n'], nelem, nnode, ndof, npdof)
fprintf(1, ['  nnodf = %9i  ndepn = %9i  nincr = %9i', ...
            '  noutp = %9i\n'] , nnodf, ndepn1, nincr, noutp)
fprintf(logfile, ['\n  nelem = %9i  nnode = %9i  ndof  = %9i', ...
            '  npdof = %9i\n'], nelem, nnode, ndof, npdof);
fprintf(logfile, ['  nnodf = %9i  ndepn = %9i  nincr = %9i', ...
            '  noutp = %9i\n'] , nnodf, ndepn1, nincr, noutp);
%
if nreac
  close all
end
%
%
% Incremental procedure
% ---------------------
%
% initialize arrays
u = zeros(ndof, 1);
du = zeros(ndof,1);
fext = spalloc(ndof, 1, 100);
fint = zeros(ndof, 1);
%
%
% initialize deformation history
hist = phist([], elem, nelem, eltp, geom, matr, data);
%
dbase = [datf,'knl'];
dbcheck = exist([dbase,'.mat']);
%
% write model data to output file
if (any(outi == 0)) && (noutp > 0) && (dbcheck == 0)
  goutp(elem, nelem, node, nnode, eltp, geom, matr, data, dof, ndof, ...
        hist, u, du, outp, noutp, outi, nouti, datf, 0, 0, 0, logfile)
end
%
iiter = 0;
iincr = 0;
conv = Inf;
dconv = 1;
refine = 0;
incrloop = 1;
niter = 0;
citer = 2;
nreduce = 0;

ntransO = 0;
cx_p = 0;

% volume average strain and stress
VASTSN = [];
%
% Get previous results if restart
%
if restart
  if dbcheck == 2
    eval(['load ',dbase]);
    fprintf(1,'\n  Loading previous mtfem-kernel \n');
    fintc = fint;
    histc = hist;
  end
end
fext_reac = zeros(ndof, 1);

% increment loop
while (incrloop)
  
  iincr = iincr + 1;
  
  aa = find(phase_vec <= iincr);
  ntransN = size(aa,2);
  
  if ntransN > ntransO
      
      % allocate nodal degrees of freedom
      [pdof, npdof, dof_n, ndof, D, depn, ndepn, cx1] = gdof(elem, nelem, eltp, nnode, pdof1, npdof1, depn1, ndepn1, ...
                                                             iincr, phase_vec, ntransO, cx_p);
      
      if ntransN == 1
          cx_p = 0;
      elseif ntransN > ntransO
          cx_p = cx1;
      end
      
      for ielem = 1:nelem
          % element type
          eltpe = eltp(elem(ielem,1),:);
          % global nod and dof numbers
          nnodee = feval(eltpe, 'nnodee');
          ndofne = feval(eltpe, 'ndofne');
          inode = elem(ielem, 4+(1:nnodee));
          
          dofe = dof(inode, :)';
          maxndofn = size(dofe, 1);
          idofn = (1:maxndofn)';
          if all(all(ndofne<=1))
              idof = dofe(ndofne==1)';
          else
              idof = dofe(idofn(:, ones(1, nnodee)) <= ...
                                            ndofne(ones(1, maxndofn), :))';
          end
       
          dofe_n = dof_n(inode, :)';
          maxndofn_n = size(dofe_n, 1);
          idofn_n = (1:maxndofn_n)';
          if all(all(ndofne<=1))
              idof1 = dofe_n(ndofne==1)';
          else
              idof1 = dofe_n(idofn_n(:, ones(1, nnodee)) <= ...
                                          ndofne(ones(1, maxndofn_n), :))';
          end
          
          u_n((idof1)) = u((idof));
          
          fint_n((idof1)) = fint((idof));
      end
      
      u = u_n';
      up = u_n';
      dof = dof_n;
      fint = fint_n';
  else
      up = u;
  end
%
% load factor and change of load factor
  incri = incr(iincr);
  ti = time(iincr); % physical time
  if iincr == 1
    Dincri = incri;
    dt = ti;
  else
    Dincri = incri - incr(iincr-1);
    dt = abs(ti - time(iincr-1)); % time is linked to the incrementation, absolute value to avoid negative times
  end
  
% temperature factor
  [incriT, DincriT] = TemEff(elem, nelem, node, eltp, incrTt, iincr, nodes_t);
  
  if incrTt(iincr) ~= 0 && dt == 0
      error('dt = 0 !!!\n');
  end
  
% show increment number and load factor
  fprintf(1, '\n  iincr = %9i  incri = %11.6g\n', iincr, incri)
  fprintf(logfile, '\n  iincr = %9i incri = %11.6g\n', iincr, incri);
%
% adjust force vectors
  fext = gfext_1(elem, nelem, node, eltp, geom, matr, ...
                              data, dof, ndof,  hist, u, du, dt, incriT, DincriT, Dincri, nodf,nnode); % thermal external force vector
%   fext = gfext(nodf, nnodf, incri, nnode, dof, ndof); % standard external force vector
  
  % get viscoplastic force increment Borg Eqs. 27 + 28 last terms on RHS
  fvp = gfvp(elem, nelem, node, eltp, geom, matr, ...
                              data, dof, ndof,  hist, u, dt);

  % residual force increment from increment of external and increment of vp
  % force
  fres = fext + fvp;
  if any(isnan(fres))
    eval(['save ',dbase ,' u fint fext_reac hist reacti iincr sxx' ]);
    error('NaN issues in fres !!!!!\n')
  end
%   fres = gfext(nodf, nnodf, Dincri, nnode, dof, ndof) + fvp;
%
% iteration loop
  iiter = 0;
  conv = Inf;
  convold = Inf;
  dconv = 1;
  dconvchk = 0;
  while ((conv > tolr) && (iiter < maxniter) && (dconv > 0))
%
%   show iteration number
    iiter = iiter + 1;
    fprintf(1, '  iiter = %9i', iiter)
    fprintf(logfile, '  iiter = %9i \n', iiter);
%
%   generate tangential stiffness matrix
    load_upd = 'y'; % does not seem to be used
    K = gk(elem, nelem, node, eltp, geom, matr, data, dof, ndof, hist, u, du);
%   store K to compute reaction forces
    Kreac = K;
% 
%   eliminate dependent degrees of freedom
    if ndepn > 0
      idofd = 1:ndepn;
      idofi = (ndepn+1):ndof;
      K = K(idofi, idofi) + D' * K(idofd, idofi) ...
          + K(idofi, idofd) * D + D' * K(idofd, idofd) * D;
      fres = fres(idofi, :) + D' * fres(idofd, :);
    end
%
%   partition for prescribed dofs
    if npdof > 0
      idofp = 1:npdof;
      idoff = (npdof+1):(ndof-ndepn);
      fres = fres(idoff, :);
      if (iiter == 1)
        fres = fres - K(idoff, idofp) * Dincri * pdof(:, 3);
      end
      K = K(idoff, idoff);
    end
%      K,
%      fres,
%   solve system

    du = K \ fres;
    
    %size(find(abs(K)<1e-100))
    %corr = find(abs(du)<1e-27);
    %du(corr) = 0;
    if any(isnan(du))
        % check if the system is solved correctly or not
        eval(['save ',dbase ,' u fint fext_reac hist reacti iincr sxx' ]);
	error('NaN issues again !!!!!\n')
	%man_du = find(isnan(du));
	%du(man_du) = zeros(size(man_du));
    end

%         trunc = 22;
%         [U,S,V] = svd(K);
%         v=diag(S);
%         Sinv=diag(1./v);
%         Vt=V(:,1:trunc);
%         Ut=U(:,1:trunc);
%         invK=Vt*Sinv*Ut';
%      du = invK * fres   

%   add prescribed dofs
    if (npdof > 0)
      if (iiter == 1)
        du = [Dincri * pdof(:, 3); du];
      else
        du = [zeros(npdof, 1); du];
      end
    end
%
%   add dependent dofs
    if ndepn > 0
      du = [D * du; du];
    end

%   construct new solution
    up = up + du;
%
%   compute reaction forces
%
    % delta fvp is put on the LHS to keep the structural reaction force
    % only
    fext_reac = fext_reac + Kreac*du-fvp;
%
%   calculate new internal nodal forces
    load_upd = 'y';
    % compute stresses and update plastic strain quantities
    [fint, hist] = gfint(elem, nelem, node, eltp, geom, matr, data, ...
                         dof, ndof, hist, up, du, dt, incriT, DincriT);
%
%   residual load vector
 %   fres = fext - fint;
    nfres = size(fext,1);
    fres = zeros(nfres,1);
%
%   determine and show convergence parameter
    conv = gconv(fres, fext, ndof, D, ndepn); % trick to always converge
    fprintf(1, '  conv  = %10.5g\n', conv)
    if iiter > citer
      dconv = convold-conv;
    end
    convold = conv;
    if ~(conv < Inf)
      dconv = -1 ;
    end
  end
%
%
  if ((iiter <= maxniter) && (dconv > 0) && (conv < tolr))
    niter = iiter;
    u = up;
%   process deformation history
    hist = phist(hist, elem, nelem, eltp, geom, matr, data);
    fintc = fint;
    histc = hist;    
    %
    if (iincr/10 - fix(iincr/10)) == 0
       % save kernel at every 10 increments for restarting
       eval(['save ',dbase ,' u fint fext_reac hist reacti iincr sxx' ]);
    end
%
% Output reaction forces / displacements

    if ndepn > 0
      n_fext = size(fext,1);
      fext_post = fext_reac;
      idofd = 1:ndepn;
      idofi = (ndepn+1):ndof;
      fext_post(idofi) = fext_reac(idofi) + D' * fext_reac(idofd);
    else
      fext_post = fext_reac;  
    end

    if nreac == 1
       disp = u(dof(reac(1,1),reac(1,2)))';
       forc = fext_post(dof(reac(1,1),reac(1,2)))';
       %   to measure the true stress - for the benchmark in Borg
       disp2 = u(dof(reac(1,1),2))';
       ssxx = fext_post(dof(reac(1,1),reac(1,2)))'/(1+disp2);
       stxx = log(1+u(dof(reac(1,1),reac(1,2)))');
       if nreac > 1
          for ir = 2:nreac
            disp = [disp u(dof(reac(ir,1),reac(ir,2)))'];
            forc = [incri * nodf(:, 3)];
          end
       end
       if ~exist('reacti')
           reacti=[disp forc];
           sxx = [stxx ssxx];
       else
           reacti=[reacti ; disp forc];
           sxx = [ sxx ; stxx ssxx];
       end
       %reacf(reacti)
       %drawnow
    else  
      reacti = [];
    end
%     
    % Output reaction forces / displacements
%
%     if nreac
%       disp = u(dof(reac(1,1),reac(1,2)))';
%       forc = fint(dof(reac(1,1),reac(1,2)))';
%       if nreac > 1
%         for ir = 2:nreac
%           disp = [disp u(dof(reac(ir,1),reac(ir,2)))'];
%           forc = [forc fint(dof(reac(ir,1),reac(ir,2)))'];
%         end
%       end
%       if ~exist('reacti')
%           reacti=[disp forc];
%       else
%         reacti=[reacti ; disp forc];
%       end
%       reacf(reacti)
%       drawnow
%     else  
%       reacti = [];
%     end

%

    % output volume average stress/strain (VASTSN)

    % get volume average stress and strain from the elements
    [vstress,vstrain,vol] = gvolav(elem, nelem, node, eltp, geom, matr, data, dof, ndof, hist, u, du);
    VASTSN = [ VASTSN ; vstress',vstrain',vol'];

%   generate output 
    if (any(iincr == outi)) && (noutp > 0)
      goutp(elem, nelem, node, nnode, eltp, geom, matr, data, dof, ...
          ndof, hist, u, du, outp, noutp, outi, nouti, datf, iincr, ...
          nreac, reacti, logfile)
    end
%     if restart
%       eval(['save ',dbase ,' u fint hist reacti nincr incr iincr outi nouti' ]);
%     end      
  end
%  
  if (dconv <= 0) || (iiter >= maxniter && conv > tolr)
    fprintf(1, '\n  Reducing load step\n')
    nreduce = nreduce + 1;
    iincr = iincr - 1;
    incr=[incr(1:iincr),(incr(iincr)+incr(iincr+1))/2,incr(iincr+1:nincr)];
    nincr=length(incr);
    nouti=length(outi);
    noutii=nouti-nnz(outi>=iincr);
    outi(noutii+1:nouti)=outi(noutii+1:nouti)+1;
    fint = fintc;
    hist = histc;
    refine = refine+1;
    if refine == 20
      fprintf(1,'\n  Excessive refining of the load step ! ')
      error('  No convergence, use other solution-algorithm !\n')
    end
  else 
    refine = 0;
  end
%
  incrloop = (iincr < length(incr));
  
  ntransO = ntransN;
  
%   if iincr == iincr_profile
%     
%     p = profile('info');
%     profsave(p,'profile_results')
%     profile off
%   end
end
%
%  write output reaction forces / displacements to file
if nreac
  fprintf(1,'\n writing reaction forces versus displacements to file. \n')
  eval(['save ', datf, 'FU.mat reacti ;'])
  
  fprintf(1,'\n writing true stress versus true strain to file. \n')
  eval(['save ', datf, 'sxx.mat sxx ;'])
end
%     p = profile('info');
%     profsave(p,'profile_results')

  fprintf(1,'\n writing  volume average stress/strain to file. \n')
  eval(['save ', datf, 'VASTSN.mat VASTSN ;'])

%
%
fprintf(1, '\n\n\n')
fprintf(logfile, '\n\n\n');
fclose(logfile);
%
% reacf(reacti)
%
%
% figure(2) 
% vm = [];
% ep = [];
% for i = [0:nincr] 
%     if exist(sprintf('%s%03d.mat', 'mTEST1', i)) ~= 2 
%         %error(sprintf('Data file `%s%04d.mat'' not found.', datf, i)) 
%         break
%     end 
%     eval(sprintf('load %s%03d', 'mTEST1', i)) 
%     vm = [vm; vonmis(2)*1e6];
%     ep = [ep; plshst(2)];
% end
% hold on
% plot (ep,vm/360,'o-','LineWidth',1.5)
% %legend ('Hardening Law','Numerical Model');
% xlabel('\epsilon^{pl}','FontSize',14,'FontName','Times')
% ylabel('\sigma_{vM} / \sigma_0','FontSize',14,'FontName','Times')
% 
% grid on

