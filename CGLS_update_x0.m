function [x,n_CGLS] = CGLS_update_x0(b,A,x0,maxit)

% The function computes the approximate CGLS update

[m,n]   = size(A);
b       = b - A*x0;
w       = zeros(n,1);  % Initialize
d       = b - A*w;
r       = A'*d;
p       = r;
aux0    = r'*r;
j       = 0;

maxit   = m;
G_prev  = inf;
discr   = sqrt(m); 

if norm(d) < discr
    iterate = false;
else
    iterate = true;
end
while (j<=maxit) && iterate
      j     = j + 1;
      d_old = d;
      y     = A*p;
      alpha = aux0/(y'*y);
      w_old = w;
      w     = w + alpha*p;       % New update
      d     = d - alpha*y;
      % Checking the discrepancy condition for the large least squares
      G     = norm(A*w - b)^2 + norm(w)^2;
      GG(j) = G;
      if norm(d) < discr || norm(d-d_old)<1e-6 || G>G_prev
          if j>1
              % If j=1, the first iteration is accepted
              w = w_old;
              j = j-1;
          end
         % The discrepancy limit is reached; stop iteration
         iterate = false;
         %  disp(['Discrepancy limit attained at iteration j =' num2str(j)]);
     end
     G_prev = G;
     r      = A'*d;
     aux1   = r'*r;   
     beta   = aux1/aux0;
     p      = r + beta*p; 
     aux0   = aux1;
end
n_CGLS = j;
x      = x0 + w; % Update with approximate IAS