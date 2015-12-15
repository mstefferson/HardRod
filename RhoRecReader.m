%% RhoRecReader.m %%
Nx   = 64;
Ny   = 64;
Nm   = 64;
Nrec = 11;

filename = 'DiffOutFT.txt';
DiffImp = importdata(filename);
%DiffImp = importdata(filename);


% fidi = fopen( filename );
% data = textscan(fidi, '%s%s', 'Delimiter',',');
DiffRec = zeros(Nx,Ny,Nm,11);

for t = 1:Nrec
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nm
                
                DiffRec(i,j,k,t) = DiffImp( j + (i-1) * Ny + (t-1) * Nx * Ny,k);
                
            end
        end
    end
end


%%
xHolder = 1;
yHolder = 1;
mHolder = Nm / 2 + 1;
tHolder   = Nrec;

Vec = reshape( DiffRec(xHolder,yHolder,:,tHolder), [1,Nm] )
figure()
plot( 1:Nm, Vec )
