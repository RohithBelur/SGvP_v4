function [incriT, DincriT] = TemEff(elem, nelem, node, eltp, incrTt, iincr, nodes_t)
  % temperature factor
  incriT = zeros(nelem,1);
  DincriT = zeros(nelem,1);
  
  if size(incrTt,1) == 1 && nodes_t(1,1) == 0
      incriT = zeros(nelem,1);
      DincriT = zeros(nelem,1);
  else
      for ielem = 1:nelem
        %   element type
        eltpe  = eltp(elem(ielem, 1), :);
        %
        %   global node and dof numbers
        nnodee = feval(eltpe, 'nnodee');
        inode  = elem(ielem, 4+(1:nnodee));
        %
        for it = 1:size(nodes_t,2)
            if ismember(inode,nodes_t(:,it)) == [1 1 1 1 1 1]
                incriT(ielem,1) = incrTt(it,iincr);
                if iincr == 1
                    DincriT(ielem,1) = incriT(ielem,1);
                else
                    DincriT(ielem,1) = incriT(ielem,1) - incrTt(it,iincr - 1);
                end
%                 it
    %         else
    %             incriT(ielem,1) = incrTt(size(incrTt,1),iincr);
    %             if iincr == 1
    %                 DincriT(ielem,1) = incriT(ielem,1);
    %             else
    %                 DincriT(ielem,1) = incriT(ielem,1) - incrTt(size(incrTt,1) ,iincr - 1);
    %             end
            end
        end
      end
  end
end

