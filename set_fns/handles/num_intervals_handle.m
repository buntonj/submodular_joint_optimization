function handle = num_intervals_handle(param)
handle = @(I) ( num_intervals(I,param) );
end