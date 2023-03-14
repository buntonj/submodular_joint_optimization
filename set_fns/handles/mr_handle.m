function handle = mr_handle(param)
handle = @(I)( modified_range(I,param) );
end