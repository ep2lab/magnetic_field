%{
Compendium of all tests for the function/class in the name of this file.
You can run the tests by executing runtests. You must add the package to
your path first. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170314
%----------------------------------------------------------------------
%}
function tests = uniform_2d_test
    tests = functiontests(localfunctions);
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_field_creation(t)
    f = magnetic_field.uniform_2d;
    f = magnetic_field.uniform_2d('B0',2);
    f = magnetic_field.uniform_2d('is_axi',0);
end
function test_field_method_calling(t)
    f = magnetic_field.uniform_2d;
    f.set_B0;
    f.set_B0(2);
    [psi,Bz,Br] = f.field_2d(0,0);
    f.psi_2d(1,1);
    Bz = f.Bz_2d(1,1);
    Br = f.Br_2d(1,1);
    B = f.B_2d(1,1);
    assert(B==sqrt(Bz^2 + Br^2));
    [psi,Bz,Br] = f.field_2d(rand(3,4),rand(3,4));    
    kB = f.curvature_2d(1,1);
    B_vector = f.B_vector_2d(1,1);
    b_vector = f.b_vector_2d(1,1);
    n_vector = f.n_vector_2d(1,1);
    assert(norm(B_vector/norm(B_vector) - b_vector) < 1e-14)
    [Z,R] = meshgrid(linspace(0,10),linspace(0,10));
    [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = f.derivatives_2d(Z,R);
    [dBz_dz,dBz_dr,dBr_dz,dBr_dr] = f.derivatives_2d(1,0);
end 
function test_plotting_method_calling(t)
    f = magnetic_field.uniform_2d;
    hf = figure(1);
    f.plot_generator_2d('color','r','linewidth',4);   
    f.plot_2d('var','Bz');
    f.plot_2d('var','Bz-Br','linestyle','none');
    close(hf);
end  
function test_streamline_propagation(t)   
    f = magnetic_field.uniform_2d;
    [Z,R] = f.streamline_2d(rand(3)*2,rand(3)*2,0.5,20);  
    [Z,R] = f.next_point_2d(rand(3)*2,rand(3)*2,0.5);
end
