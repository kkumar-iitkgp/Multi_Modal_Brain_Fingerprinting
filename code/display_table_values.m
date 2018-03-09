function display_table_values(data_matrix, cell_column, cell_row_label)

sTable = array2table(data_matrix,'RowNames',cell_row_label,'VariableNames',cell_column);

sTable

end