% This is the source code for the experiment of comparing the solution time of 
% kOCP described in the paper Higher Order Cone Program. 
% We try to solve the problem 
% minimize   tr(X)
% subject to X \in K_k^n
%            tr(A_iX)=f_i, i =1,2,...,p
% where K_k^n is the kth order cone. 
% We choose n =15 and p =5 and A_i, f_i are randomly generated. 
% We compare the solution time for different k. 
% We run several trials and average the solution time for each k.
% In each trial, we run k from 1 to 15 and make the A_i, f_i fixed.
% The way we solve it is using the following equivalent program 
% minimize   tr(\sum_{i_1<\dots <i_k} M^{i_1\dots i_k})
% subject to tr(A_i\sum_{i_1<\dots <i_k} M^{i_1\dots i_k})=f_i ,i=1,...,5
%            T_{i_1\dots i_k} ?M^{i_1\dotsi_k} positive semidefinite
% where T_{i_1\dots i_k} is the operator that take only the k*k submatrix
% with index (i,j),i,j in \{i_1,...i_k\},see more description in theorem
% 3.3 of the paper.
% The solution time is the number after "Interior-point optimizer
% terminated. Time:" and is recorded by hand. 

cvx_solver MOSEK  % The default solver is SDPT3 which is slow when k is large
tt=3 % total number of different trials. 
n=10;%; %dimension of matrix variable X
M = zeros(n,tt); % matrix record the optimal value.
for  q = 1:tt
q  
p=5; % number of constraints
v= 1:n; % integer 1 to n
AA = zeros(n,n,p); % constraint matrix

for i = 1:p
    AA(:,:,i) = rand(n);
    AA(:,:,i) = triu(AA(:,:,i))+ triu(AA(:,:,i),1);
end
% the above is for building the constraint amtrix

a=rand(p,1);% constraint value
B=1:b;

k=1 %1OCP or Diagonally Dominant Matrix Program (DDP)
cvx_begin 
      variable X(n,n) symmetric 
      minimize trace(X)
      subject to
          for e =1:n
              X(e,e)>=sum(abs(X(e,[1:e-1 e+1:end])));
          end  
          for P=1:p
          trace(AA(:,:,P)*X) == a(P);
          end  
cvx_end
M(1,q)= cvx_optval;

for k = 2:n-1 %kOCP
k
C = nchoosek(v,k);% all choose k combinations of 1:n ;
b = nchoosek(n,k);%combination number n choose k;
cvx_begin sdp
      variable m(k,k,b) symmetric semidefinite
      expression Y(n,n)
      for i = 1:b
          Y(C(i,:),C(i,:)) = Y(C(i,:),C(i,:)) + m(:,:,i);
      end 
      minimize trace(Y)
      subject to
          for P=1:p
          trace(AA(:,:,P)*Y) == a(P);
          end
cvx_end
M(k,q)= cvx_optval;
end 

n
cvx_begin sdp %nOCP or SDP 
      variable X(n,n) symmetric semidefinite
      minimize trace(X)
      subject to
          for P=1:p
          trace(AA(:,:,P)*X) == a(P);
          end 
cvx_end
M(n,q)= cvx_optval;
end 
