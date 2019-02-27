%% Compile Mex files
%

% Copyright (c) 2017-2019, SeHyoun Ahn
% License: 2-clause BSD license. You can check the license at
%   <https://sehyoun.com/sparse_grid/license.html> XXXX

mex -v -largeArrayDims COMPFLAGS='-O2' find_neighbor_mex.c;
mex -v -largeArrayDims COMPFLAGS='-O2' union_sorted.c;
