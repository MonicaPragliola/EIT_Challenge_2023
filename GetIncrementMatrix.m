function L_full = GetIncrementMatrix(CondVertices,CondTopol,Edges);

% This program computes the finite difference matrix. We assume that the
% boundary values of the perturbation in conductivity are zero.
%--------------------------------------------------------------------------
% Calls to: none
% D Calvetti, E Somersalo
%--------------------------------------------------------------------------

% Increment matrix, computing increments over edges. L_full contains
% boundary edges that are later removed

n_vert = size(CondVertices,2);
n_edge = size(Edges,1);

rows   = [1;1]*(1:n_edge); % Two non-zero entries in each row
rows   = rows(:);
cols   = Edges';               % Non-zero enries in columns pointing to the end verices
cols   = cols(:);
vals   = [1;-1]*ones(1,n_edge); % values +1 and -1
vals   = vals(:);
L_full = sparse(rows,cols,vals,n_edge,n_vert);