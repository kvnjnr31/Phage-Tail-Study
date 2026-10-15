%Apply reflecting B.C.'s ----------------------------------
    if abs(z_new) >= z_upper
   
        z_new = z_old-(z_new-z_old);                % z position.
        z_old = z_upper*z_new./abs(z_new);
    else
        z_old = z_new;
    end