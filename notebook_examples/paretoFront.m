function [params,front] = paretoFront(Xs,Fs,eFs)
% Return the Pareto frontier of a set of points
% 
% ----------
% Parameters
% ----------
% Xs: N x k matrix
%     N samples with the k parameters of each sample
% Fs: N x m matrix
%     N samples with the m values of each objective function
% eFs: N x m matrix
%     N samples with the error in each objective function
% 
% -------
% Returns
% -------
% params: p x k matrix
%         parameters of the p Pareto optimal points
% front: p x m matrix
%         values of the p Pareto optimal points

N = size(Xs,1);
params = [];
front = [];
for i = 1:N
   if paretoPoint(Fs,i)
       params = [params; Xs(i,:)];
       front = [front; Fs(i,:)];
   end
end