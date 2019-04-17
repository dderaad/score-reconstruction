function filter = Shannonfilter(t,a,tau)
%Shannonfilter
%   Shannon-type 'harsh' filter (takes values of either 0 or 1)
filter = abs(t-tau) <= a/2;

end

