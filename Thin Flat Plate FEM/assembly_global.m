
function global_matrix = assembly_global(local_matrix, n)
M = local_matrix;
n_total = 3 * n^2;
e = n - 1;
dofs_element = length(local_matrix);
M_g = zeros(n_total);



for j=1:e
    for i=1:e
        indexes = get_indexes(i,j,n,dofs_element);
        M_g(indexes,indexes) = M_g(indexes,indexes) + M;
    end
end


global_matrix = M_g;

