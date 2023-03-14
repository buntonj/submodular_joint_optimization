% visualize the fit we made with DD regression

i = 2; % pick the product number to visualize

[x, I] = sort(avg_price(i,:),'ascend');
%x = avg_price;
max_price = x(end);
y = num_sales;
y_fit = zeros(num_intervals,num_products);
for j = 1:num_intervals
    y_fit(j,:) = avg_price(:,j)'*beta + alpha;
end

scatter(x,y(i,:));
hold on;
scatter(x,y_fit(:,i));