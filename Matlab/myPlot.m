A=csvread('../out.ppq.allpatients.v20191109.csv',1,0);

num_patients = size(A,1)/672;
num_hours=size(A,1)/num_patients;
xx = A(1:num_hours,2);
C=[];

for pid=0:(num_patients-1)
    
    B = A( A(:,1)==pid, : );
    C=[C, B(:,3)];
    %plot( B(:,2), B(:,3), 'k-' ); hold on; 
    %plot( B(:,2), B(:,4), 'b-' ); 
end

qq=quantile(C, [.025 .25 .5 .75 .975], 2);

subplot(1,2,1)
plot( xx, log10( qq(:,3) ), 'b-', 'LineWidth', 2 ); hold on;
plot( xx, log10( qq(:,2) ), 'b-' ); 
plot( xx, log10( qq(:,4) ), 'b-' ); 
plot( xx, log10( qq(:,1) ), 'b--' ); 
plot( xx, log10( qq(:,5) ), 'b--' ); 

grid on;
xlabel('HOURS');
ylabel('log10( mg PPQ )');
axis([0 336 0 3.5]);


%%

A=csvread('../out.ppq.allpatients.v20191110.csv',1,0);

num_patients = size(A,1)/672;
num_hours=size(A,1)/num_patients;
xx = A(1:num_hours,2);
C=[];

for pid=0:(num_patients-1)
    
    B = A( A(:,1)==pid, : );
    C=[C, B(:,3)];
    %plot( B(:,2), B(:,3), 'k-' ); hold on; 
    %plot( B(:,2), B(:,4), 'b-' ); 
end

qq=quantile(C, [.025 .25 .5 .75 .975], 2);

subplot(1,2,2)
plot( xx, log10( qq(:,3) ), 'b-', 'LineWidth', 2 ); hold on;
plot( xx, log10( qq(:,2) ), 'b-' ); 
plot( xx, log10( qq(:,4) ), 'b-' ); 
plot( xx, log10( qq(:,1) ), 'b--' ); 
plot( xx, log10( qq(:,5) ), 'b--' ); 

grid on;
xlabel('HOURS');
ylabel('log10( mg PPQ )');
axis([0 336 0 3.5]);