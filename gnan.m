function nan = gnan(m, n)
%
% function nan = gnan(m, n)
%
% Generation of m x n matrix of NaNs.
%
% input:
%   m        :  number of rows
%   n        :  number of columns
%
% output:
%   nan      :  NaN-matrix
%

%
% construct matrix
nan = NaN;
nan = nan(ones(1, m), ones(1, n));
%
%
%
