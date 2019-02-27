function [output] = union_sorted(v1, v2)
% INTERNAL FUNCTION: Fail-safe file for union_sorted.c
%
% LICENSE:
%   Copyright 2017-2019 SeHyoun Ahn
%   BSD 2-clause see <https://github.com/sehyoun/>
%

  output = sort([v1(:);v2(:)]);
  output = output([diff(output) ~= 0;true]);
end
