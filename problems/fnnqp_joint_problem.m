%JOINT SUBMODULAR-CONVEX QP PROBLEM CLASS
% special class for problems of the form:
%
% Minimize  0.5*x'*H(x)*x + b*x + F(supp(x))
% x in R^n
%
% where supp(x) denotes the set of nonzero entries in x.
% The algorithm is guaranteed to solve the problem exactly when $H$ is
% submodular, i.e. has all nonpositive off-diagonal entries.
%
% Uses the Fujishige-Wolfe MNP algorithm to find the minimum of a submodular
% function.
%
% To declare, requires H, b, and lambda parameters.
% Sample declaration and usage:
% 
% H = eye(3);
% b = -ones(3,1);
% lambda = 0.5;
% problem = joint_problem(H, b, lambda);
% I = h.optimize()
%%

classdef fnnqp_joint_problem < handle
    properties
        H
        b
        n
        V
        lb
        ub
        lambda
        TOL = 1e-9;
        F0 = 0;
        opt = sfo_opt({'verbosity_level',0});
        submod_function
        bottom
    end
    
    methods
        function self = fnnqp_joint_problem(H,b,lambda,G)%wt,function_choice)
            if any((H > 0) - diag(diag(H>0)),'all')
                %fprintf('WARNING: ')
                %fprintf('Quadratic form is not submodular. \n')
            end
            if any(eig(H) < 0)
                fprintf('WARNING: ')
                fprintf('Quadratic form is not convex. \n')
            end
            self.H = H;
            self.b = b;
            self.lambda = lambda;
            self.n = length(b);
            self.lb = zeros(self.n,1);
            self.ub = Inf(self.n,1);
            self.V = 1:self.n;
            self.bottom = [];
            self.submod_function = G;
            self.F0 = self.joint_value([]); % normalization factor
        end
        
        %function to compute value of a subset
        function val = joint_value(self,I)
            %set the correct support constraint
            Aeq = eye(self.n);
            for i = 1:length(self.bottom)
                Aeq(self.bottom(i),self.bottom(i)) = 0;
            end
            for i = 1:length(I)
                Aeq(I(i),I(i)) = 0;
            end
            [~,fval] = fnnqp(self.H(I,I),self.b(I),1e-10);
            val = fval + self.lambda*self.submod_function(I) - self.F0;
        end
        
        function x_star = set_to_minimizer(self,I)
            x_star = zeros(self.n,1);
            [x_star(I),~] = fnnqp(self.H(I,I),self.b(I),1e-10);
        end
        
        % use semi-gradient methods from Iyer, Jegelka, Bilmes to prune the
        % lattice [0, V] into the sublattice [bottom, top]
        % with a new ground set V = top \ bottom
        % runs in O(n) time and offers a ~50% search space reduction.
        function prune_space(self)
            [self.bottom, ~, self.V] = semigradient_pruning(@(I) self.joint_value(I),self.V);
        end
        
        %perform optimization problem
        function [I_star,x_star,subopt,costs] = optimize(self)
            [I_star,subopt,costs] = bootleg_min_norm_point(sfo_fn_wrapper(@(A) self.joint_value(A)),self.V,self.opt);
            I_star = unique([self.bottom,I_star]);%find(abs(x_star) > self.TOL);
            x_star = self.set_to_minimizer(I_star);
        end
    end
    
    
    
    methods (Static)
        %default is cardinality
        function F = weighted_cardinality(I,params)
            F = sum(params.w(I));
        end
        
        function F = sqrt_cardinality(I,params)
            F = sqrt(sum(params.w(I)));
        end
        
        function F = modified_range(I,params)
            % Modified range function + cardinality
            if any(I)
                F = length(params.w) + max(I) - min(I) + length(sfo_unique_fast(I));
            else
                F = 0;
            end
        end
        
        function F = num_intervals(I,params)
            % number of sets of contiguous nonzero entries
            I = sort(I);
            d = (I-circshift(I,-1)) ~= -1;
            F = sum(d(1:(length(d)-1))) + sum(params.w(I)) + 1;
        end
        
        % grouped penalties
        % requires params to have:
        % G : n x num_groups binary array
        % w : num_groups x 1 weighting array
        function F = groups(I,params)
            F = sum(params.w(any(params.G(I,:),1)));
        end
        
        %{
        function F = modified_range(I,w)
            % weighted cardinality
            %F = sum(w(I));
            
            % pure cardinality
            %F = length(sfo_unique_fast(I));
            
            % Modified range function + cardinality
            %if any(I)
            %    F = length(w) + max(I) - min(I) + length(sfo_unique_fast(I));
            %else
            %    F = 0;
            %end
            
            % number of sets of contiguous nonzero entries
            I = sort(I);
            d = (I-circshift(I,-1)) ~= -1;
            F = sum(d(1:(length(d)-1))) + sum(w(I)) + 1;
        end
        %}
    end
end