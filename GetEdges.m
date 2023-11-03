function [E] = GetEdges(T)

%Return the edges of the triangulation (number of edges x 2)
% Input: Topology matrix T, n_triangle x 3 integer array
% Output: E List of edges, a n_edges x 2 integer array
%--------------------------------------------------------------------------
% D Calvetti, E Somersalo
% Calls to: None
%--------------------------------------------------------------------------

dim=size(T);
nedge=3;
e1=T(1,1:2);
e2=T(1,2:3);
e3=[T(1,1), T(1,3)];
E=[e1;e2;e3];
       
for i=2:dim(1)
        e1=T(i,1:2);
        e2=T(i,2:3);
        e3=[T(i,1), T(i,3)];
        not_new1=(any(find((E(:,1)==e1(1)) & (E(:,2)==e1(2)))) || any(find((E(:,1)==e1(2)) & (E(:,2)==e1(1)))));
        not_new2=(any(find((E(:,1)==e2(1)) & (E(:,2)==e2(2)))) || any(find((E(:,1)==e2(2)) & (E(:,2)==e2(1)))));
        not_new3=(any(find((E(:,1)==e3(1)) & (E(:,2)==e3(2)))) || any(find((E(:,1)==e3(2)) & (E(:,2)==e3(1)))));
        if (~not_new1)
            nedge=nedge+1;
            E=[E;e1];
        end
        if (~not_new2)
            nedge=nedge+1;
            E=[E;e2];
        end
        if (~not_new3)
            nedge=nedge+1;
            E=[E;e3];
        end
end
end

