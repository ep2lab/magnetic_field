%{
Compendium of all tests for the function/class in the name of this file.
You can run the tests by executing runtests. You must add the package to
your path first. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170314
%----------------------------------------------------------------------
%}
function tests = loop_2d_test
    tests = functiontests(localfunctions);
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_loop_creation(t)
    loop = magnetic_field.loop_2d;
    loop = magnetic_field.loop_2d('RL',23);
    loop = magnetic_field.loop_2d('ZL',23);
    loop = magnetic_field.loop_2d('I',23);        
end
function test_field_method_calling(t)
    loop = magnetic_field.loop_2d;
    loop.set_B0;
    loop.set_B0(2);
    [psi,Bz,Br] = loop.field_2d(0,0);
    loop.psi_2d(1,1);
    Bz = loop.Bz_2d(1,1);
    Br = loop.Br_2d(1,1);
    B = loop.B_2d(1,1);
    assert(B==sqrt(Bz^2 + Br^2));
    [psi,Bz,Br] = loop.field_2d(rand(3,4),rand(3,4));    
    kB = loop.curvature_2d(1,1);
    B_vector = loop.B_vector_2d(1,1);
    b_vector = loop.b_vector_2d(1,1);
    n_vector = loop.n_vector_2d(1,1);
    assert(all(B_vector/norm(B_vector) == b_vector))
    [Z,R] = meshgrid(linspace(0,10),linspace(0,10));
    [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = loop.derivatives_2d(Z,R);
    [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = loop.derivatives_2d(1,0);
end
function test_fastloop(t)
    loop = magnetic_field.loop_2d('interpKE',1);
    [psi,Bz,Br] = loop.field_2d(0,0);
    [psi,Bz,Br] = loop.field_2d(rand(3,4),rand(3,4));   
    [Z,R] = meshgrid(linspace(0,10),linspace(0,10));
    [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = loop.derivatives_2d(Z,R);
    [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = loop.derivatives_2d(1,0);
end
function test_plotting_method_calling(t)
    loop = magnetic_field.loop_2d;
    hf = figure(1);
    loop.plot_generator_2d('color','r','linewidth',4);   
    loop.plot_2d('var','Bz');
    loop.plot_2d('var','Bz-Br','linestyle','none');
    close(hf);
end  
function test_streamline_propagation(t)   
    loop = magnetic_field.loop_2d;
    [Z,R] = loop.streamline_2d(rand(3)*2,rand(3)*2,0.5,20);  
    [Z,R] = loop.next_point_2d(rand(3)*2,rand(3)*2,0.5);
end
