A=csvread('../out.ppq.allpatients.v20191109.csv',1,0);

for pid=0:99
    
    B = A( A(:,1)==pid, : );
    plot( B(:,2), B(:,3), 'k-' ); hold on; 
    plot( B(:,2), B(:,4), 'b-' ); 
end


