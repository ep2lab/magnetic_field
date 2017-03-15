%{
Compendium of all tests for the function/class in the name of this file.
You can run the tests by executing runtests. You must add the package to
your path first. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170314
%----------------------------------------------------------------------
%}
function tests = loop_3d_test
    tests = functiontests(localfunctions);
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_loop_creation(t)
    loop = magnetic_field.loop_3d;
    loop = magnetic_field.loop_3d('RL',23);
    loop = magnetic_field.loop_3d('ZL',23);
    loop = magnetic_field.loop_3d('I',23);        
end
function test_field_method_calling(t)
    loop = magnetic_field.loop_3d;
    loop.set_B0;
    loop.set_B0(2);
    [psi,Bz,Br] = loop.field_2d(0,0);
    [Bx,By,Bz] = loop.field_3d(0,0,0);
    Bx = loop.Bx_3d(0,0,0);
    By = loop.By_3d(0,0,0);
    Bz = loop.Bz_3d(0,0,0);
    B = loop.B_3d(0,0,0);
    loop.axis = [1,1,1]; % It should be normalized internally
    loop.origin = [0,0,2];
    [Bx,By,Bz] = loop.field_3d(0,0,0);    
    [Bx,By,Bz] = loop.field_3d(rand(3,4),rand(3,4),rand(3,4));    
    B_vector = loop.B_vector_3d(1,1,1);
    b_vector = loop.b_vector_3d(1,1,1);
    assert(norm(B_vector/norm(B_vector) - b_vector) < 1e-14)
    [X,Y,Z] = meshgrid(linspace(0,10),linspace(0,20),linspace(0,30));
    [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = loop.derivatives_3d(X,Y,Z);    
end  
function test_plotting_method_calling(t)
    loop = magnetic_field.loop_3d;
    hf = figure(1);
    loop.plot_generator_3d('color','r','linewidth',4);   
    loop.plot_3d('var','Bz');
    loop.plot_3d('var','Bz-Bx','linestyle','none');
    loop.plot_3d('o',[1,1,1],'var','Bz-Bx','linestyle','none');
    close(hf);
end  
function test_streamline_propagation(t)   
    loop = magnetic_field.loop_3d;
    [x,y,z] = loop.streamline_3d(rand(3)*2,rand(3)*2,rand(3)*2,0.5,20);  
    [x,y,z] = loop.next_point_3d(rand(3)*2,rand(3)*2,rand(3)*2,0.5);
end
       
