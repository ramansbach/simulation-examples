function f = paretoPointErr(Xs,eXs,i)
%Function that determines whether a point is on the Pareto frontier
%
%----------
%Parameters
%----------
% Xs: N x m matrix of objective functions for all points
% eXs: N x m matrix of error in the objective functions for all points
% i: int
%   row where point to be tested is found
%
%--------
% Returns
%--------
% f: bool
% true if xstar is Pareto optimal, false otherwise
% -----
% Notes
% -----
% A point x* is Pareto Optimal if there is no such X that F_i(x) <= F_i(x*)
% for all i = 1,...,N, that is, there is no point that is better in all
% objective functions (it's okay if it's the same point -- strict equality
% over all objectives is fine)
f = false;
N = size(Xs,1);
m = size(Xs,2);
for j = 1:N
   if j ~= i
       if dominatesPoint(Xs(j,:),Xs(i,:),eXs(i,:),eXs(j,:))
            2;
            return
       end
   end
end
f = true;
end

