%This is the function used to generate the matrix used in the first part
%of the Crank-Nicholson scheme. The matrix is a tridiagonal system with
%zero flux boundary conditions included.

%The main loop for creating the tridiagonal matrix. Each r(i) is different
%and must be accounted for - the loop takes this into consideration. Cell
%arrays are used as a simple way of producing block matrices. The B.C.s
%are also included here, in this circumstance they are zero flux.
function X=CNmatrix1(a,theta)
b=cell(1,length(a{3}));
for pp=1:length(a{3})
    [b{pp}]=gallery('tridiag',length(a{1}),(1-theta)*a{6}(pp),1-2*(1-theta)*a{6}(pp),(1-theta)*a{6}(pp));
    %Here we insert the boundary conditions which will carry over for each
    %cell.
    cellarray=[b{pp}];
    cellarray(1,1)=1-2*(1-theta)*a{6}(pp); 
    cellarray(1,2)=2*(1-theta)*a{6}(pp);
    cellarray(length(a{1}),length(a{1})-1)=2*(1-theta)*a{6}(pp); 
    cellarray(length(a{1}),length(a{1}))=1-2*(1-theta)*a{6}(pp);
    [b{pp}]=cellarray;
end
%This step concatenates the cells into a block diagonal matrix.
X=(blkdiag(b{:}));
end


