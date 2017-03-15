%{
Compendium of all tests for the function/class in the name of this file.
You can run the tests by executing runtests. You must add the package to
your path first. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170314
%----------------------------------------------------------------------
%}
function tests = array_3d_test
    tests = functiontests(localfunctions);
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_array_creation(t)
    a = magnetic_field.array_3d;
    f1 = magnetic_field.loop_3d;
    f2 = magnetic_field.loop_3d('ZL',3);
    a.generators = {f1,f2};
end
function test_field_method_calling(t) 
    f1 = magnetic_field.loop_3d;
    f2 = magnetic_field.loop_3d('ZL',3);
    f = magnetic_field.array_3d('generators',{f1,f2});
    f.set_B0;
    f.set_B0(2);
    [Bx,By,Bz] = f.field_3d(0,0,0); 
 	[Bx,By,Bz] = f.field_3d(rand(3,4),rand(3,4),rand(3,4));    
    B_vector = f.B_vector_3d(1,1,1);
    b_vector = f.b_vector_3d(1,1,1);
    assert(all(B_vector/norm(B_vector) == b_vector))
    [X,Y,Z] = meshgrid(linspace(0,10),linspace(0,20),linspace(0,30));
    [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = f.derivatives_3d(X,Y,Z);    
end 
function test_plotting_method_calling(t)
    f1 = magnetic_field.loop_3d;
    f2 = magnetic_field.loop_3d('ZL',5,'RL',6);
    f = magnetic_field.array_3d('generators',{f1,f2});
    hf = figure(1);
    f.plot_generator_3d('color','r','linewidth',4);   
    f.plot_3d('var','Bz');
    f.plot_3d('var','Bz-By','linestyle','none');
    close(hf);
end  
function test_streamline_propagation(t)   
    f1 = magnetic_field.loop_3d;
    f2 = magnetic_field.loop_3d('ZL',3);
    f = magnetic_field.array_3d('generators',{f1,f2});
    [x,y,z] = f.streamline_3d(rand(3)*2,rand(3)*2,rand(3)*2,0.5,20);  
    [x,y,z] = f.next_point_3d(rand(3)*2,rand(3)*2,rand(3)*2,0.5);
end 
