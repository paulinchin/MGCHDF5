  % Store outputs
    clear Tevol;
    clear OHevol;
    clear OIevol;
    %Set OH Model Coefficients   
    b=[0,0,0,0,0,0,0.08,0.17,0.27,0.48];
    a8=[3.9,10.5,25,25,25,25,25,25,25,25];
    a9=[ 2 0 0 0 0 0 0 0 0
         0 4 0 0 0 0 0 0 0 
         0 1 7 0 0 0 0 0 0 
         0 1 2 10 0 0 0 0 0
         0 1 2 6 16 0 0 0 0
         1 1 3 6 11 22 0 0 0
         4 6 9 12 16 23 32 0 0
         4 6 8 10 14 19 25 33 0
         28 29 31 32 34 36 38 40 42];
    a10=[0,0.58,1.0,1.7,3.0,5.2,9.1,16,70,48];
    aa=[ 22.74,0,0,0,0,0
         30.43,15.42,0,0,0,0
         28.12,40.33,2.032,0,0,0
         20.30,69.77,7.191,0.299,0,0
         11.05,99.42,15.88,1.315,0.051,0
         4,125.6,27.94,3.479,0.274,0.010
         2.34,145.1,42.91,7.165,0.847,0.063
         8.6,154.3,59.98,12.68,2.007,0.230
         23.72,148.9,78.64,19.94,4.053,0.620];
     
        A62=1.75;  A83=.78;   A31=27.9637;
        
    %Set OI Model Coefficients
    beta=0.04;
    k2=2.5e-12;
    K10=9e-14;
    zz=[7,12,16,29,33,37,38,41,45,59,88,167,274,426,647,921,1269];
    k3O2=K10.*zz;
    k3O=5.*K10.*zz;