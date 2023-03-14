function [bottom, top, V_new] = semigradient_pruning(F,V)

% GROW STEP
X = [];
while 1
    candidates = setdiff(V,X); % initialize pool of possible elements to add
    FX = F(X); % function evaluated at current set    
    X_new = X; % start with current set of elements
    for j = 1:length(candidates) % loop through possible additions
        if (F([X,candidates(j)]) - FX < 0) % if we stand to decrease F by adding j
            X_new = [X_new, candidates(j)]; % then do it
        end
    end
    
    if ~any(setdiff(X_new,X)) % if we didn't grow
        break % terminate
    end
    X = X_new; % otherwise, take a step forward
end
bottom = X_new;

X = V;
X_new = [];
% SHRINK STEP
while 1
    FX = F(X); % function evaluated at current set
    
    X_new = X; % we will remove elements from current X
    for j = 1:length(X)
        if (FX - F(X(X~=X(j)))) > 0 % if adding j to X increased F
            X_new = setdiff(X_new,X(j)); % then remove
        end
    end
    
    if ~(any(setdiff(X,X_new))) % if we didn't shrink the set this step
        break % we terminated
    end
    
    X = X_new; % we start again at the newly shrunk set
end
top = X_new;

% NEW GROUND SET
% the new lattice to search over is [bottom, top],
% which is all possible supersets of bottom and subsets of top.
V_new = setdiff(top,bottom);
end