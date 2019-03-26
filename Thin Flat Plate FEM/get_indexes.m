function out = get_indexes(i,j,n, dofs_element)
    out = zeros(1,dofs_element);
    out(1) = dof_n(i,j,n);
    out(2) = out(1) + 1 ;
    out(3) = out(2) + 1 ;
    
    out(4) = dof_n(i+1,j,n);
    out(5) = out(4) + 1 ;
    out(6) = out(5) + 1 ;
    
    out(7) = dof_n(i+1,j+1,n);
    out(8) = out(7) + 1 ;
    out(9) = out(8) + 1 ;
    
    out(10) = dof_n(i,j+1,n);
    out(11) = out(10) + 1 ;
    out(12) = out(11) + 1 ;