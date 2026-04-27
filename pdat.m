function [node, nnode, elem, nelem, eltp, geom, matr, data, ...
          pdof1, npdof1, nodf, nnodf, depn1, ndepn1, incr, nincr, ...
	      tolr, maxniter, reac, nreac, outp, noutp, outi, ...
	      nouti, restart, time, phase_vec, incrTt, nodes_t,kersave,ch_pdof,ch_vec,sv_pdof,fact] = pdat(datf);
%function [node, nnode, elem, nelem, eltp, geom, matr, data, ...
%          pdof, npdof, nodf, nnodf, depn, ndepn, incr, nincr, ...
%	   tolr, maxniter, reac, nreac, outp, noutp, outi, ...
%          nouti, restart] = pdat(datf);
%
% Reading and processing of data file.
%
% input:
%   datf     :  name of data file
%
% output:
%   node     :  nodal coordinates
%   nnode    :  number of nodes
%   elem     :  structure connectivity
%   nelem    :  number of elements in structure
%   eltp     :  element types
%   geom     :  geometry groups
%   matr     :  material property groups
%   data     :  element data groups
%   pdof     :  prescribed degrees of freedom
%   npdof    :  number of prescribed degrees of freedom
%   nodf     :  nodal forces
%   nnodf    :  number of nodal forces
%   depn     :  dependencies between degrees of freedom
%   ndepn    :  number of dependencies
%   incr     :  increment load factors
%   nincr    :  number of increments
%   tolr     :  convergence tolerance
%   maxniter :  maximum number of iterations
%   reac     :  desired output reaction forces
%   nreac    :  number of output reaction forces
%   outp     :  output items
%   noutp    :  number of output items
%   outi     :  increments for which output is to be generated
%   nouti    :  number of output increments
%   restart  :  restart-option
%   time     :  time domain


%
% Read model data
% ---------------
%
% check existence of data file
if exist(datf) ~= 2
  error(['Data file not found: ', datf, '.'])
end
%
% read data from data file
eval(datf)
%
%
% Process model data
% ------------------
%
% add absent items to model
if ~exist('geom')
  geom = [ 0 ];
  fprintf(1, '  Warning: No geometry data supplied.\n')
elseif isempty(geom)
  geom = [ 0 ];
end
if ~exist('matr')
  matr = [ 0 ];
  fprintf(1, '  Warning: No material data supplied.\n')
elseif isempty(matr)
  matr = [ 0 ];
end
if ~exist('data')
  data = [ 1 ];
  fprintf(1, '  Warning: No element data supplied.\n')
elseif isempty(data)
  data = [ 1 ];
end
if ~exist('nodf')
  nodf = [];
end
if ~exist('nodf')
  nodf = [];
end
if ~exist('depn1')
  depn1 = [];
end
if ~exist('outp')
  outp = [];
end
if ~exist('outi')
  outi = [ 0 1];
end
if ~exist('incr')
  incr = [0 1];
end
if ~exist('reac')
  reac = [];
end
if ~exist('restart')
  restart = 0 ;
end
if ~exist('tolr')
  tolr = 1.e-7 ;
  fprintf(1, '  Warning: No convergence tolerance specified, default value 1.e-7 taken.\n')  
end
if ~exist('maxniter')
  maxniter = 50 ;
  fprintf(1, '  Warning: No maximum iterations specified, default value 50 taken.\n')  
end
if ~exist('time')
  time = [];
  fprintf(1, '  Warning: No time steps specified.\n')
end
%
%  Reorder numbering sequence in elements if necessary
%
% for ielem=1:size(elem,1)
%   zer=sum(elem(ielem,5:size(elem,2))==0);
%   nodex=[elem(ielem,5:size(elem,2)-zer) elem(ielem,5)];
%   dx=node(nodex(2:length(nodex)),1)-node(nodex(1:length(nodex)-1),1);
%   ym=(node(nodex(1:length(nodex)-1),2)+node(nodex(2:length(nodex)),2))/2.;
%   vol=sum(dx.*ym);
%   if vol>0
%     elem(ielem,6:size(elem,2)-zer)=fliplr(elem(ielem,6:size(elem,2)-zer));
%   end
% end
% 
% 
% % condens out double specified rows in pdof
% cnd = [0];
% keep = ones(size(pdof,1),1);
% for i=1:size(pdof,1)
%   if ~any(cnd==i)
%     a=pdof(i,1); b=pdof(i,2); c=pdof(i,3);
%     exc = find((a==pdof(:,1)).*(b==pdof(:,2)).*(c==pdof(:,3)));
%     cnd = [cnd; exc];
%     for j=1:length(exc)
%       if exc(j)~=i
%         keep(exc(j)) = 0; 
%       end 	
%     end  
%   end  
% end  
% pdof = pdof(logical(keep),:);
% %
% % condens out double specified rows in ndof
% cnd = [0];
% keep = ones(size(nodf,1),1);
% for i=1:size(nodf,1)
%   if ~any(cnd==i)
%     a=nodf(i,1); b=nodf(i,2); c=nodf(i,3);
%     exc = find((a==nodf(:,1)).*(b==nodf(:,2)).*(c==nodf(:,3)));
%     cnd = [cnd; exc];
%     for j=1:length(exc)
%       if exc(j)~=i
%         keep(exc(j)) = 0; 
%       end 	
%     end  
%   end  
% end  
% nodf = nodf(logical(keep),:);

%
% determine length of a number of arrays
nnode = size(node, 1);
nelem = size(elem, 1);
npdof1 = size(pdof1, 1);
nnodf = size(nodf, 1);
ndepn1 = size(depn1, 1);
nincr = length(incr);
noutp = size(outp, 1);
nouti = length(outi);
nreac = size(reac,1);
%
% display warnings
if noutp == 0
  fprintf(1, '  Warning: No output items specified.\n')
end
%
%
%
