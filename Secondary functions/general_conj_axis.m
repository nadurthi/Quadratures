function D=general_conj_axis(n,m)
%n is the dimension of system
%m is the number of ones 1(s) I need in each vector
dr=prod_conjugate_dir(m);
C = nchoosek(1:n,m);
[rc,cc]=size(C);
[rdr,cdr]=size(dr);
% D=zeros(2^m*nchoosek(n,m),n);
% for i=0:rc:rdr*rc-rc
%     for j=1:1:rc
%         for k=1:1:m
%         D(i+j,C(j,k))=dr(floor(i/rc)+1,k);
% %         D(i+j,C(j,2))=dr(floor(i/rc)+1,2);
% %         D(i+j,C(j,3))=dr(floor(i/rc)+1,3);
%         end
%     end
% end
D=[];
for i=1:1:nchoosek(n,m)
    A=zeros(2^m,n);
    for j=1:1:m
        A(:,C(i,j))=dr(:,j);
    end
    D=vertcat(D,A);
end