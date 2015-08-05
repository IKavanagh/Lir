function [Z, N, M] = zmread(filename)
% Needs expanded comment
%zmread(filename) Reads a complex matrix from a file pointed to by filename
% where the first line of the file is the rank of the matrix and the subsequent
% rank lines are the dimensions of the matrix.

fid = fopen(filename);
rank = fscanf(fid, '%d\n', 1);
N = fscanf(fid, '%d ', rank); fscanf(fid, '\n');

if rank == 1
    N(2) = 1;
    rank = 2;
end

N = reshape(N, 1, rank);

fprintf(1, 'Reading matrix of rank %d with dimensions ', rank);
for i = 1:rank-1
    fprintf(1, '%dx', N(i));
end
fprintf(1, '%d\n', N(end));

Z = zeros(N);
M = prod(N);

for i = 1:M
    Z(i) = fscanf(fid, '%f', 1) + 1i*fscanf(fid, '%fi\n', 1);
end

fclose(fid);

end