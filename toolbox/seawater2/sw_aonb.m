
function [AONB]= aonb(S, PTMP, P)

% SW_AONB    Calculate alpha/beta (a on b)
%================================================================
% SW_AONB   $Revision: 1.4 $   $Date: 1994/10/10 04:31:38 $
%           Copyright (C) CSIRO, Nathan Bindoff 1993
%
% USAGE: [AONB]= aonb(S, PTMP, P)
%
% DESCRIPTION
%    Calculate alpha/beta.  See sw_alpha.m and sw_beta.m
%
% INPUT:  (all must have same dimensions)
%   S    = salinity              [psu      (PSS-78) ]
%   PTMP = potential temperature [degree C (IPTS-68)]
%   P    = pressure              [db]
%          (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT
%   AONB  = alpha/beta [psu/degree_C]
%
% AUTHOR:   N.L. Bindoff  1993
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE: 
%    McDougall, T.J. 1987. "Neutral Surfaces"
%    Journal of Physical Oceanography vol 17 pages 1950-1964, 
%
% CHECK VALUE: 
%    aonb=0.34763 psu C^-1 at S=40.0 psu, ptmp=10.0 C, p=4000 db
%================================================================

% Modifications
% 93-04-22. Phil Morgan,  Help display modified to suit library
% 93-04-23. Phil Morgan,  Input argument checking

% CHECK INPUT ARGUMENTS
if nargin~=3
  error('sw_aonb.m: requires 3 input arguments')
end %if

T=PTMP;

% CHECK S,T,P dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
[mp,np] = size(P);

  
% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   error('check_stp: S & T must have same dimensions')
end %if

% CHECK OPTIONAL SHAPES FOR P
if     mp==1  & np==1      % P is a scalar.  Fill to size of S
   P = P(1)*ones(ms,ns);
elseif np==ns & mp==1      % P is row vector with same cols as S
   P = P( ones(1,ms), : ); %   Copy down each column.
elseif mp==ms & np==1      % P is column vector
   P = P( :, ones(1,ns) ); %   Copy across each row
elseif mp==ms & np==ns     % PR is a matrix size(S)
   % shape ok 
else
   error('check_stp: P has wrong dimensions')
end %if
[mp,np] = size(P);
 

  
% IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
Transpose = 0;
if mp == 1  % row vector
   P       =  P(:);
   T       =  T(:);
   S       =  S(:);   

   Transpose = 1;
end %if
%***********************************************************
%End $RCSfile: check_stp.Minc,v $
%***********************************************************

% BEGIN
	 c1=fliplr([ 0.665157e-1, 0.170907e-1, ...
	    -0.203814e-3, 0.298357e-5, ...
            -0.255019e-7]);
         c2=fliplr([ 0.378110e-2, ...
            -0.846960e-4]);
         c2a=fliplr([0.0 -0.164759e-6, ...
            -0.251520e-11]);
         c3=[-0.678662e-5];
         c4=fliplr([+0.380374e-4, -0.933746e-6, ...
            +0.791325e-8]);
         c5=[0.512857e-12];
         c6=[-0.302285e-13];
%
% Now calaculate the thermal expansion saline contraction ratio adb
%
        [m,n] = size(S);
        sm35  = S-35.0*ones(m,n);
        AONB  = polyval(c1,PTMP) + sm35.*(polyval(c2,PTMP)...
               + polyval(c2a,P)) ...
               + sm35.^2*c3 + P.*polyval(c4,PTMP) ...
               + c5*(P.^2).*(PTMP.^2) + c6*P.^3;

return
%----------------------------------------------------------------------

