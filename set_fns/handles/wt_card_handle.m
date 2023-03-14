function handle = wt_card_handle(param)
handle = @(I) (weighted_cardinality(I,param));
end