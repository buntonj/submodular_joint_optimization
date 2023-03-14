%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        PROJECTED GRADIENT DESCENT METHOD               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_star,x_star,hist] = pgm_solve(prob,opt_params)
maxiter = opt_params.maxiter;
if isfield(opt_params,'thresh')
    thresh = opt_params.thresh;
else
    thresh = 1e-8;
end
if isfield(opt_params,'step')
    step = opt_params.step;
else
    step = 5e-2;
end
n = prob.n;
I_vect = (rand(n,1));
I_set = find(I_vect);
last = prob.joint_value(I_set);
hist(1) = last;
for i = 1:maxiter
    subgrad = sfo_polyhedrongreedy(sfo_fn_wrapper(@(A) prob.joint_value(A)),1:n,I_vect);
    I_vect = I_vect - step*subgrad';
    I_vect(I_vect > 1) = 1;
    I_vect(I_vect < 0) = 0;
    I_set = find(abs(I_vect)>eps);
    new = prob.joint_value(I_set);
    hist(i+1) = new;
    %if new-last <= eps
    %    break
    %end
    last = new;
end
I_star = I_set;
Aeq = eye(n);
for i = 1:length(I_star)
    Aeq(I_star(i),I_star(i)) = 0;
end
[x_star,~] = cplexqp(prob.H,prob.b,[],[],Aeq,zeros(prob.n,1));
I_star = find(x_star);
end