%This is a function for calculating the kinetics to be inserted into the
%main calculation loop. I take each individual term within the kinetics and
%calculate these seperately, and then concatenate them at the end. This is
%the function for universal mutation.

function B=kinetics(a,dt)

%Initialise the arrays to be used within the loop.
b=cell(3,1);
c=cell(1,1);
d=cell(length(a{2}));
%This cell is for storing the mutation matrix with 0 on the diagonals for
%calculating the final terms within the kinetics easier.
c{2}=a{5}-diag(diag(a{5}));
%This loop calculates the kinetics for each individual species.
for ii=1:length(a{2})
    %This term is the competition coefficient part.
    b{1}=(-1)*sum(bsxfun(@times,a{4}(ii,:),a{1}),2);
    %This term is the positive incoming mutations.
    b{2}=sum(bsxfun(@times,c{2}(:,ii)',a{1}),2);
    %This term is the negative leaving mutations.
    b{3}=(-1)*sum(c{2}(ii,:))*a{1}(:,ii);
    %This term is the full kinetic equation using the above.
    d{ii}=(a{2}(ii).*a{1}(:,ii).*(1+b{1})+b{2}+b{3})*dt;
end
%Concatenating the matrix into the form used for solving numerically.
B=cat(1,d{:});

end



