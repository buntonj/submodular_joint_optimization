function handle = groups_handle(param)
handle = @(I) (groups(I,param));
end