if exist('num_products')~= 1 % check if the problem dimension is already set
    num_products = 10; % number of products to consider
    fprintf('No input, running with n = %d \n',num_products);
end
products = unique(data{:,3}); % pulls the names of first num_products
indices = randsample(seed,size(products,1),num_products,false);
%indices = 1:num_products;
products = products(indices);
clear indices;
%indices = randsample(size(data,1),num_products,false); % randomly sample indices for products
%products = unique(data{indices,3});

num_intervals = 3; % number of sales intervals to consider be careful, because no sales in an interval causes errors

f = false(size(data,1),1);
for i = 1:num_products
    f = or(f, strcmp(products{i},data.Description));
end

filtered_data = data(f,:);
start_date = min(filtered_data.InvoiceDate);
end_date = max(filtered_data.InvoiceDate);

interval_length = (end_date - start_date)/num_intervals; % divide time into even intervals

num_sales = zeros(num_products,num_intervals);
avg_price = zeros(num_products,num_intervals);

for i=1:num_products
    f1 = strcmp(products{i},filtered_data.Description);
    for j = 1:num_intervals
        int_start = start_date + (j-1)*interval_length;
        int_end = start_date + j*interval_length;
        f2 = (filtered_data.InvoiceDate >= int_start) & (filtered_data.InvoiceDate <= int_end);
        num_sales(i,j) = sum(filtered_data(f1 & f2,:).Quantity);
        avg_price(i,j) = mean(filtered_data(f1 & f2,:).UnitPrice);
        if isnan(avg_price(i,j))
            avg_price(i,j) = 0;
        end
    end
end

% fix zeros in average price
for i = 1:num_products
    for j = 1:num_intervals
        if avg_price(i,j) == 0
            % search backward
            for k = j-1:-1:1
                if avg_price(i,k) > 0
                    avg_price(i,j) = avg_price(i,k);
                    break
                end
            end
            
            % if you found a price to pull, exit
            if avg_price(i,j) > 0
                continue
            end
            
            % otherwise, search forward for a price to pull
            for k = j:num_intervals
                if avg_price(i,k) > 0
                    avg_price(i,j) = avg_price(i,k);
                    break
                end
            end
        end
    end
end

% taking out the larger trash
clear f;
clear f1;
clear f2;
clear filtered_data;
%clear data; % you're done with that bigass data file, get rid of it