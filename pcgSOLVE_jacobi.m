function [Q] = pcgSOLVE_jacobi(K,F, CG_tol_power)
% add to the line     du =pcgSOLVE( K , fres); ~ line 165
% TOLERANCE
%NUMDDL=size(F,1);
%if NUMDDL<10^6
%    CHOPRE=1e-1;
%    PCGPRE=1e-6;
%else
%    CHOPRE=0.5e-1;
%    PCGPRE=1e-6;
%end

% FACTORIZATION
% try

%alpha = max(sum(abs(K),2)./diag(K))-2;
%L=ichol(K, struct('type','ict','droptol',CHOPRE,'diagcomp',alpha));
%[Q flag relres iter resvec]=pcg(K,F,PCGPRE,10000,L,L');
%disp(' ');
%disp(['PCG TIME : ' num2str(toc) ' s']);
%save('resvec','resvec');

[Q, flag, relres, iter]=pcg(K,F,(10^(CG_tol_power)),10000,diag(diag(K)));
%save('resvec','resvec');

if flag~=0
    disp(' ');
    disp('!!! PCG NON-ZERO FLAG');
    if     flag==1 disp('!!! Maximum iteration reached');
    elseif flag==2 disp('!!! Matrix ill-conditionned');
    elseif flag==3 disp('!!! Pcg stagnated');
    elseif flag==4 disp('!!! Too large or too small quantity');
    end
end
%iter
%relres
%resvec
