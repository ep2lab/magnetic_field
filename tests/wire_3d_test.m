%{
Compendium of all tests for the function/class in the name of this file.
You can run the tests by executing runtests. You must add the package to
your path first. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170314
%----------------------------------------------------------------------
%}
function tests = wire_3d_test
    tests = functiontests(localfunctions);
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_wire_creation(t)
    f = magnetic_field.wire_3d; 
    f = magnetic_field.wire_3d('I',23);    
end
function test_field_method_calling(t)
    f = magnetic_field.wire_3d;
    f.set_B0;
    f.set_B0(2);
    [psi,Bz,Br] = f.field_2d(0,0);
    [Bx,By,Bz] = f.field_3d(0,0,0);
    Bx = f.Bx_3d(0,0,0);
    By = f.By_3d(0,0,0);
    Bz = f.Bz_3d(0,0,0);
    B = f.B_3d(0,0,0);
    f.direction = [0,1,0]; % It should be normalized internally
    f.point = [0,0,0]; 
    B = f.B_2d([1,2,3]-2,[1,2,3]-3.5); % For manual comparison
    B_vector = f.B_vector_3d(1,1,1);
    b_vector = f.b_vector_3d(1,1,1);
    assert(norm(B_vector/norm(B_vector) - b_vector) < 1e-14)
    [X,Y,Z] = meshgrid(linspace(0,10),linspace(0,20),linspace(0,30));
    [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = f.derivatives_3d(X,Y,Z);    
end  
function test_advanced_comparison_2d_3d(t)
    f = magnetic_field.wire_3d;
    f.set_B0;  
    f.direction = [0,1,0]; % It should be normalized internally
    f.point = [0,0,0]; 
    x = rand(20,40);
    z = rand(20,40);
    [Bx,By,Bz] = f.field_3d(x,x*0,z);    
    [psi,Bz2,Bx2] = f.field_2d(z,x);  
    max(abs(Bx(:)-Bx2(:)))
    max(abs(Bz(:)-Bz2(:)))
end
function test_plotting_method_calling(t)
    f = magnetic_field.wire_3d;
    hf = figure(1);
    f.plot_generator_3d('color','r','linewidth',4);   
    f.plot_3d('var','Bz');
    f.plot_3d('var','Bz-Bx','linestyle','none');
    f.plot_3d('o',[1,1,1],'var','Bz-Bx','linestyle','none');
    close(hf);
end  
function test_streamline_propagation(t)   
    f = magnetic_field.wire_3d;
    [x,y,z] = f.streamline_3d(rand(3)*2,rand(3)*2,rand(3)*2,0.5,20);  
    [x,y,z] = f.next_point_3d(rand(3)*2,rand(3)*2,rand(3)*2,0.5);
end
       
