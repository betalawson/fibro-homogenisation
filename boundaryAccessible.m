function Vfrac = boundaryAccessible(occ)
% This function takes as input a matrix of logical values, 1's
% corresponding to obstruction and 0's to conductive material (pore space),
% then determines what proportion of the material is accessible from the
% boundaries. This gives a "true" volume fraction, that ignores
% inaccessible pores.

% Read out the size of the input matrix
[Ny, Nx] = size(occ);

% Create a list of boundary locations (linear indexing)
bounds_i = []; bounds_j = [];
% Left boundary
bounds_i = [bounds_i, ones(1,Ny)];
bounds_j = [bounds_j, 1:Ny];
% Right boundary
bounds_i = [bounds_i, Nx*ones(1,Ny)];
bounds_j = [bounds_j, 1:Ny];
% Top boundary
bounds_i = [bounds_i, 2:Nx-1];
bounds_j = [bounds_j, ones(1,Nx-2)];
% Bottom boundary
bounds_i = [bounds_i, 2:Nx-1];
bounds_j = [bounds_j, Ny*ones(1,Nx-2)];

bounds = [];
bounds = [bounds, 1 : 1 : Ny ];                     % Left boundary
bounds = [bounds, ((Nx-1)*Ny+1) : 1 : (Ny*Nx) ];    % Right boundary  
bounds = [bounds, 2 : Ny : ((Nx-2)*Ny+1) ];         % Top boundary
bounds = [bounds, Ny*2 : Ny : (Nx-1)*Ny ];          % Bottom boundary 

% Initialise a matrix that will be used to track what is accessible
A = double(occ);

% Apply a method that converts all accessed 0's in the input matrix into
% 2's. Final result will be composed of 0's for unaccessed material, 1's
% for obstacles, 2's for accessed material
for k = 1:length(bounds)
   
    % Initialise the method at this boundary site (if the boundary site has
    % already been reached in the course of processing the algorithm
    % initialised at a separate boundary site, the algorithm will simply
    % terminate immediately without continuing)
    A = exploreAccessible( A, bounds_i(k), bounds_j(k) );
    
end

% Finally, calculate the Vfrac, as the number of accessible sites divided
% by the total
Vfrac = sum(A(:) == 2) / numel(A);


end



% function A = exploreAccessible(A, i, j)
% % This recursive function repeatedly moves to accessible neighbours not yet
% % labelled as accessible, and labels them accessible.
% 
% % Only process if the current site has not already been accessed
% if A(j,i) == 0
%     
%     % Convert this site to accessible, then process all of its neighbours
%     A(j,i) = 2;
%     
%     if i > 1  % Can check left
%         
%             A = exploreAccessible(A, i-1, j);
%         
%         if j > 1 % Can check left-down
%                 A = exploreAccessible(A, i-1, j-1);
%         end
%         if j < Ny % Can check left-up
%                 A = exploreAccessible(A, i-1, j+1);
%         end
%         
%     end
%     
%     if i < Nx  % Can check right
%         
%             A = exploreAccessible(A, i+1, j);
%         
%         if j > 1 % Can check right-down
%                 A = exploreAccessible(A, i+1, j-1);
%         end
%         if j < Ny % Can check right-up
%                 A = exploreAccessible(A, i+1, j+1);
%         end
%         
%     end
%     
%     if j > 1 % Can check down
%             A = exploreAccessible(A, i, j-1);
%     end
%     
%     if j < Ny % Can check up
%             A = exploreAccessible(A, i, j+1);
%     end
%     
% end
% 
% end
% 

function A = exploreAccessible(A, i, j)
% This recursive function repeatedly moves to accessible neighbours not yet
% labelled as accessible, and labels them accessible.

% Only process if the current site has not already been accessed
if A(j,i) == 0
    
    % Convert this site to accessible, then process all of its neighbours
    A(j,i) = 2;
    
    % Couch all attempts to explore neighbours in try statements, so that
    % there is no error thrown when trying to look outisde the domain
    try
        A = exploreAccessible(A, i-1, j);
    end
    try
        A = exploreAccessible(A, i+1, j);
    end
    try
        A = exploreAccessible(A, i-1, j-1);
    end
    try
        A = exploreAccessible(A, i-1, j+1);
    end
    try
        A = exploreAccessible(A, i, j-1);
    end
    try
        A = exploreAccessible(A, i, j+1);
    end
    try
        A = exploreAccessible(A, i+1, j-1);
    end
    try
        A = exploreAccessible(A, i+1, j-1);
    end
    
end

end

