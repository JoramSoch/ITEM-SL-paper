function r = MD_unirnd(a, b, N)
% _
% Random Numbers from a Continuous Uniform Distribution
% FORMAT r = MD_unirnd(a, b, N)
% 
%     a - the lower end of the uniform distribution
%     b - the upper end of the uniform distribution
%     N - the number of samples to be drawn
% 
%     r - an N x 1 vector of uniform random numbers
%
% FORMAT r = MD_unirnd(a, b, N) draws N values from the continuous uniform
% distribution with minimum a and maximum b. This is done by sampling r
% from the standard uniformation distribution x ~ U(0,1) and then
% calculating y = x*(b-a) + a [1].
% 
% References:
% [1] Wikipedia: "Uniform distribution (continuous)";
%     URL: https://en.wikipedia.org/wiki/Uniform_distribution_(continuous).
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 16/10/2018, 09:55
%  Last edit: 31/10/2023, 13:03


% Generate random numbers
%-------------------------------------------------------------------------%
r = (rand(N,1) * (b-a)) + a;