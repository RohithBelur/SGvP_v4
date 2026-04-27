mex -c -g -O -largeArrayDims findDet.F90
mex -c -g -O -largeArrayDims gdndxi.F90
mex -c -g -O -largeArrayDims gintp.F90
mex -c -g -O -largeArrayDims gn.F90
mex -c -g -O -largeArrayDims inv.F90

mex -c -g -O -largeArrayDims ghardening.F90 inv.o

% mex -c -g -O -largeArrayDims plbq2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% mex -g -O -largeArrayDims plbq2t15.F90 plbq2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o

mex -c -g -O -largeArrayDims plbr2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
mex -g -O -largeArrayDims plbr2t15.F90 plbr2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o

% mex -c -g -O -largeArrayDims plta2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% mex -g -O -largeArrayDims plta2t15.F90 plta2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o

% mex -c -g -O -largeArrayDims plte2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% mex -g -O -largeArrayDims plte2t15.F90 plte2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o

mex -c -g -O -largeArrayDims plt12t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
mex -g -O -largeArrayDims plt12t15.F90 plt12t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o

mex -c -g -O -largeArrayDims plt22t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
mex -g -O -largeArrayDims plt22t15.F90 plt22t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o

mex -c -g -O -largeArrayDims plt32t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
mex -g -O -largeArrayDims plt32t15.F90 plt32t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o

mex -c -g -O -largeArrayDims plt42t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
mex -g -O -largeArrayDims plt42t15.F90 plt42t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o

mex -c -g -O -largeArrayDims plgm2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
mex -g -O -largeArrayDims plgm2t15.F90 plgm2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% 
% mex -c -g -O -largeArrayDims plnt2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% mex -g -O -largeArrayDims plnt2t15.F90 plnt2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% 
% mex -c -g -O -largeArrayDims plnq2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% mex -g -O -largeArrayDims plnq2t15.F90 plnq2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% 
% mex -c -g -O -largeArrayDims plaa2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% mex -g -O -largeArrayDims plaa2t15.F90 plaa2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% 
% mex -c -g -O -largeArrayDims plbb2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% mex -g -O -largeArrayDims plbb2t15.F90 plbb2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o

% mex -c -g -O -largeArrayDims pltr2t15Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% mex -g -O -largeArrayDims pltr2t15.F90 pltr2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o ghardening.o
% 
% mex -c -g -O -largeArrayDims plcs2t12Fct.F90 findDet.o gdndxi.o gintp.o gn.o inv.o
% mex -g -O -largeArrayDims plcs2t12.F90 plcs2t12Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o
% 
% mex -c -g -O -largeArrayDims ingkTnrFct.F90
% mex -g -O -largeArrayDims ingkTnr.F90 ingkTnrFct.o plnn2t15Fct.o plnt2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o
% 
% mex -c -g -O -largeArrayDims gfintTnrFct.F90
% mex -g -O -largeArrayDims gfintTnr.F90 gfintTnrFct.o plnn2t15Fct.o plnt2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o
%
% mex -c -g -O -largeArrayDims find.F90
% 
% mex -c -g -O -largeArrayDims ingdofTFct.F90
% mex -g -O -largeArrayDims ingdofT.F90 ingdofTFct.o plnn2t15Fct.o plnt2t15Fct.o findDet.o gdndxi.o gintp.o gn.o inv.o find.o
