function [vstressGR, vstrainGR, vol_tGR, vstressT1, vstrainT1, vol_tT1, vstressT2, vstrainT2, vol_tT2, vstressT3, vstrainT3, vol_tT3, vstressT4, vstrainT4, vol_tT4, vstressMo, vstrainMo, vol_tMo] =  gvolav(elem, nelem, node, eltp, geom, matr, data, dof, ndof, hist, u, du)
%
% 
% Computing Volume Average of Stress/Strain through elements.
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
%  vstress   :  Volume Average of Stress
%  vstrain   :  Volume Average of Strain
%


% in all simulations plane strain modeling assumption
vstressGR=zeros(4,1);
vstrainGR=zeros(3,1);
vol_tGR=0;

vstressT1=zeros(4,1);
vstrainT1=zeros(3,1);
vol_tT1=0;

vstressT2=zeros(4,1);
vstrainT2=zeros(3,1);
vol_tT2=0;

vstressT3=zeros(4,1);
vstrainT3=zeros(3,1);
vol_tT3=0;

vstressT4=zeros(4,1);
vstrainT4=zeros(3,1);
vol_tT4=0;

vstressMo=zeros(4,1);
vstrainMo=zeros(3,1);
vol_tMo=0;

%  get element data in element loop
for ielem = 1:nelem 
  
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
    if strcmp(eltpe(1:2),'pl') && strcmp(eltpe(3:4),'cs')
        vstresse = 0;
        vstraine = 0;
        volumee  = 0;
    else
        % get element volume average stress/strain
        [vstresse,vstraine,volumee] = feval(eltpe, 'stsnGR', nodee, geome, matre, int64(datae), histe, ue, due,[0]);
        
        %  restor average of each element in the total value of the structure
        vstressGR=vstressGR+vstresse'; % element volume x stress in the element (average of the IPs)
        vstrainGR=vstrainGR+vstraine';
        vol_tGR = vol_tGR + volumee;
        
        [vstresse,vstraine,volumee] = feval(eltpe, 'stsnT1', nodee, geome, matre, int64(datae), histe, ue, due,[0]);
        
        %  restor average of each element in the total value of the structure
        vstressT1=vstressT1+vstresse'; % element volume x stress in the element (average of the IPs)
        vstrainT1=vstrainT1+vstraine';
        vol_tT1 = vol_tT1 + volumee;
        
        [vstresse,vstraine,volumee] = feval(eltpe, 'stsnT2', nodee, geome, matre, int64(datae), histe, ue, due,[0]);
        
        %  restor average of each element in the total value of the structure
        vstressT2=vstressT2+vstresse'; % element volume x stress in the element (average of the IPs)
        vstrainT2=vstrainT2+vstraine';
        vol_tT2 = vol_tT2 + volumee;
        
        [vstresse,vstraine,volumee] = feval(eltpe, 'stsnT3', nodee, geome, matre, int64(datae), histe, ue, due,[0]);
        
        %  restor average of each element in the total value of the structure
        vstressT3=vstressT3+vstresse'; % element volume x stress in the element (average of the IPs)
        vstrainT3=vstrainT3+vstraine';
        vol_tT3 = vol_tT3 + volumee;
        
        [vstresse,vstraine,volumee] = feval(eltpe, 'stsnT4', nodee, geome, matre, int64(datae), histe, ue, due,[0]);
        
        %  restor average of each element in the total value of the structure
        vstressT4=vstressT4+vstresse'; % element volume x stress in the element (average of the IPs)
        vstrainT4=vstrainT4+vstraine';
        vol_tT4 = vol_tT4 + volumee;
    
        [vstresse,vstraine,volumee] = feval(eltpe, 'stsnMo', nodee, geome, matre, int64(datae), histe, ue, due,[0]);
        
        %  restor average of each element in the total value of the structure
        vstressMo=vstressMo+vstresse'; % element volume x stress in the element (average of the IPs)
        vstrainMo=vstrainMo+vstraine';
        vol_tMo =vol_tMo + volumee;
    end
end

%get global volume average stress/strain
%vstress=vstress/vol_t;
%vstrain=vstrain/vol_t;