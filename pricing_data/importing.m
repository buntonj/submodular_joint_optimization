%filename = 'online_retail.xlsx';
filename = 'online_retail_II.xlsx';
data = readtable(filename);
%data = readtable(filename2);
%col_names = data.Properties.VariableNames;
%data2.Properties.VariableNames = col_names;
%data = vertcat(data,data2);
%clear data2;
data.Properties.VariableNames = {'InvoiceNo','StockCode','Description','Quantity','InvoiceDate','UnitPrice','CustomerID','Country'};
data(data.Quantity <= 0,:) = []; % remove rows with zero or negative sale entries indicating returns
data(data.UnitPrice <= 0,:) = []; % remove rows with zero or negative sale price indicating returns
clear filename;