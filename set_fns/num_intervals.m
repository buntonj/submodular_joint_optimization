function F = num_intervals(I,params)
    % number of sets of contiguous nonzero entries
    I = sort(I);
    d = (I-circshift(I,-1)) ~= -1;
    F = sum(d(1:(length(d)-1))) + sum(params.w(I)) + 1;
end