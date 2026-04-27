function [xi, w] = gintp(iplc, nintp)
%
% function [xi, w] = gintp(iplc, nintp)
%
% Generation of integration point information: local coordinates and
% weight factors for numerical integration.
%
% input:
%   iplc     :  interpolation code
%   nintp    :  number of integration points
%
% output:
%   xi       :  local coordinates of the integration points
%   w        :  corresponding weight factors
%

%
% One-dimensional
% ---------------
%
if iplc(1) == 'b'
%
%   1 integration point
  if nintp == 1
    xi =               [  0 ];
    w  =               [  2 ];
%
%   2 integration points
  elseif nintp == 2
    xi = (1/sqrt(3)) * [ -1
                          1 ];
    w  =               [  1
                          1 ];
%
%   3 integration points
  elseif nintp == 3
    xi = sqrt(3/5) *   [ -1
                          0
                          1 ];
    w  = (1/9) *       [  5
                          8
                          5 ];
  elseif nintp == 4
    xi =  [ -0.861136311594953
            -0.339981043584856
            +0.339981043584856
            +0.861136311594953 ];
    w  =  [ 0.347854845137454
            0.652145154862546
            0.652145154862546
            0.347854845137454 ];
%
%   illegal number of ingegration points
  else
    error(sprintf('Illegal number of integration points: %g.', ...
          nintp))
  end
%
%
% Two-dimensional, triangular
% ---------------------------
%
elseif iplc(1) == 't'
%
%   1 integration point
  if nintp == 1
    xi = (1/3) *       [  1   1   1 ];
    w  =               [  1 ];
%
%   3 integration points
  elseif nintp == 3
    xi = (1/2) *       [  1   1   0
                          0   1   1
                          1   0   1 ];
    w  = (1/3) *       [  1
                          1
                          1 ];
%
%   4 integration points
  elseif nintp == 4 
    xi = (1/15) *      [  5   5   5
                          9   3   3
                          3   9   3
                          3   3   9 ];
    w  = (1/48) *      [ -27
                          25
                          25
                          25 ];
%   6 integration points
  elseif nintp == 6
    al1 = 0.816847572980459;
    be1 = 0.091576213509771;
    al2 = 0.108103018168070;
    be2 = 0.445948490915965; 
    xi =  [  al1   be1   be1
             be1   al1   be1
             be1   be1   al1
             al2   be2   be2
             be2   al2   be2
             be2   be2   al2 ];
    w  =  [  0.109951743655322
             0.109951743655322
             0.109951743655322
             0.223381589678011
             0.223381589678011
             0.223381589678011 ];
%
%   7 integration points
  elseif nintp == 7
    al1 = 0.059715871789770;
    be1 = 0.470142064105115;
    al2 = 0.797426985353087;
    be2 = 0.101286507323456;
    xi =               [ 1/3   1/3   1/3
                         al1   be1   be1
                         be1   al1   be1
                         be1   be1   al1
                         al2   be2   be2
                         be2   al2   be2
                         be2   be2   al2 ];
    w  =               [ 0.225030000300000
                         0.132394152788506
                         0.132394152788506
                         0.132394152788506
                         0.125939180544827
                         0.125939180544827
                         0.125939180544827 ];
%
%   illegal number of integration points
  else
    error(sprintf('Illegal number of integration points: %g.', ...
          nintp))
  end
%
%
% Two-dimensional, quadrilateral
% ------------------------------
%
elseif iplc(1) == 'q'
%
%   1 integration point
  if nintp == 1
    xi =               [  0   0 ];
    w  =               [  4 ];

%
%   4 integration points (2x2)
  elseif nintp == 4
    xi = (1/sqrt(3)) * [ -1  -1
                          1  -1
                          1   1
                         -1   1 ];
    w  =               [  1
                          1
                          1
                          1 ];
%
%   9 integration points (3x3)
  elseif nintp == 9
    xi = sqrt(3/5) *   [ -1  -1
                          0  -1
                          1  -1
                          1   0
                          1   1
                          0   1
                         -1   1
                         -1   0
                          0   0 ];
    w  = (1/81) *      [  25
                          40
                          25
                          40
                          25
                          40
                          25
                          40
                          64 ];
%
%   illegal number of integration points
  else
    error(sprintf('Illegal number of integration points: %g.', ...
          nintp))
  end
%
%
% Lobatto integration scheme, quadratic
% -------------------------------------
%
elseif(iplc(1) == 'l')
%   1 integration point
	if nintp == 1
        xi = [  0 ];
        w  = [  2 ];
%
%   2 integration points
	elseif nintp == 2
        xi = [ -1 ; 1 ];
        w  = [  1 ; 1 ];
%
%   3 integration points
	elseif nintp == 3
        xi = [ -1 ; 0 ; 1 ];
        w  = (1/3) * [ 1 ; 4 ; 1 ];
	elseif nintp == 4
        xi =  [ -.447213592499958 ; -1 ; 1 ; .447213592499958 ];
        w  =  (1/6) * [ 5 ; 1 ; 1 ; 5 ];
%
%   illegal number of ingegration points
	else
        error(sprintf('Illegal number of integration points: %g.', ...
          nintp))
	end
    %  
%
%
% Illegal interpolation code
% --------------------------
%
else
  error(sprintf('Illegal interpolation code : %s.', iplc))
end
%
%
%
