%{
Compendium of all tests for the function/class in the name of this file.
You can run the tests by executing runtests. You must add the package to
your path first. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170314
%----------------------------------------------------------------------
%}
function tests = path_3d_test
    tests = functiontests(localfunctions);
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_path_creation(t)
    f = magnetic_field.path_3d; 
    f = magnetic_field.path_3d('I',23);  
    theta = linspace(0,2*pi,50);
    x = 1*cos(theta);
    y = 2*sin(theta);
    z = theta*0;
    f = magnetic_field.path_3d('I',23,'x',x,'y',y,'z',z);
end
function test_field_method_calling(t)
    theta = linspace(0,2*pi,20);
    x = 1*cos(theta);
    y = 2*sin(theta);
    z = theta*0;
    f = magnetic_field.path_3d('I',23,'x',x,'y',y,'z',z);
    f.set_B0;
    f.set_B0(2); 
    [Bx,By,Bz] = f.field_3d(0,0,0);
    Bx = f.Bx_3d(0,0,0);
    By = f.By_3d(0,0,0);
    Bz = f.Bz_3d(0,0,0);
    B = f.B_3d(0,0,0);  
    B_vector = f.B_vector_3d(1,1,1);
    b_vector = f.b_vector_3d(1,1,1);
    assert(norm(B_vector/norm(B_vector) - b_vector) < 1e-14)
    [X,Y,Z] = meshgrid(linspace(0,10),linspace(0,20),linspace(0,30));
    [dBx_dx,dBx_dy,dBx_dz,dBy_dx,dBy_dy,dBy_dz,dBz_dx,dBz_dy,dBz_dz] = f.derivatives_3d(X,Y,Z);    
end   
function test_plotting_method_calling(t)
    f = magnetic_field.path_3d;
    hf = figure(1);
    f.plot_generator_3d('color','r','linewidth',4);   
    f.plot_3d('var','Bz');
    f.plot_3d('var','Bz-Bx','linestyle','none');
    f.plot_3d('o',[1,1,1],'var','Bz-Bx','linestyle','none');
    close(hf);
end  
function test_streamline_propagation(t)   
    f = magnetic_field.path_3d;
    [x,y,z] = f.streamline_3d(rand(3)*2,rand(3)*2,rand(3)*2,0.5,20);  
    [x,y,z] = f.next_point_3d(rand(3)*2,rand(3)*2,rand(3)*2,0.5);
end
       
