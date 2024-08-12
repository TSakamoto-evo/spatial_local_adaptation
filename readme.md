Code for Sakamoto (2024)

## Folder "est_prob"
Codes for stepping-stone simulation with one patchy region to calculate establishment probability.

Each subfolder contains codes for different spatial dimension.

### Usage
1. Compile files using the following command:
```
g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out
```

2. Run the exective file by command:
```
./XXX.out arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9 arg10 arg11 arg12
```
```arg1``` (int): population density, $\rho$  
```arg2``` (int): the number of segment divisions per unit length, $n$  
```arg3``` (double): migrational variance, $\sigma^2$  
```arg4``` (double): selection coefficient in environment 1, $s_1$  
```arg5``` (double): selection coefficient in environment 2 (with inversed sign), $s_2$  
```arg6``` (double): radius of the target patch, $R$  
```arg7``` (int): the total number of replications  
```arg8``` (int): accessory number (no effect on simulation behavior)  
```arg9``` (int): accessory number (no effect on simulation behavior)  
```arg10``` (int): distance between the mutation origin and the center of the target patch (in a unit of grid length $1/n$)  
```arg11``` (int): It should be $0$ when you start new simulations. If there is a suspended data, you can specify the number of established cases in the suspended data and restart the simulation.  
```arg12``` (int): It should be $0$ when you start new simulations. If there is a suspended data, you can specify the number of replications in the suspended data and restart the simulation.  

### Output file
Final file: "est_prob.txt"  
1st col: distance between the mutation origin and the center of the target patch (in an absolute scale, not a unit of $1/n$), $r$  
2nd col: $\sigma^2$  
3rd col: $s_1$  
4th col: $s_2$  
5th col: the number of runs in which allele A is established  
6th col: total number of runs  
7th col: ```arg8```  
8th col: ```arg9```   

## Folder "wait_time"
Codes for stepping-stone simulation with two patchy regions to calculate waiting time until local adaptation through immigration.

Each subfolder contains codes for different spatial dimension.

### Usage
1. Compile files using the following command:
```
g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out
```

2. Run the exective file by command:
```
./XXX.out arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9 arg10 arg11
```
```arg1``` (int): population density, $\rho$  
```arg2``` (int): the number of segment divisions per unit length $n$  
```arg3``` (double): migrational variance, $\sigma^2$  
```arg4``` (double): selection coefficient in environment 1, $s_1$  
```arg5``` (double): selection coefficient in environment 2 (with inversed sign), $s_2$  
```arg6``` (double): radius of the target patch, $R$  
```arg7``` (double): radius of the source patch, $R_s$  
```arg8``` (int): accessory number (no effect on simulation behavior)  
```arg9``` (int): accessory number (no effect on simulation behavior)  
```arg10``` (int): distance between the source and target patch, $D$  
```arg11``` (int): accessory number (no effect on simulation behavior)  

### Output file
Final file: "wait_time.txt"  
1st col: $D$  
2nd col: $\sigma^2$  
3rd col: $s_1$  
4th col: $s_2$  
5th col: $R_s$  
6th col: an indicator to check whether the simulation successfully finished. It should be 2 when it ended successfully. It is 0 when the number of generations exceeds maximum generation, which can be specified in the main function. It should be 1 when the allele went extinct globally during the simulation (due to too small size for the source patch).  
7th col: waiting time  
8th col: ```arg8```  
9th col: ```arg9```  
10th col: ```arg11```


## Folder "theory"
Codes for numerical calculations. First, you need to calculate the establishment probability using "calculate_est.cpp". Then, you can calculate $\lambda_{\text{mig}}$ using "calculate_waiting_time.cpp". "calculate_wait_time_RC.cpp" was used to calculate the theoretical prediction presented in Ralph and Coop (2015).

### "calculate_est.cpp"
#### Usage
1. Set parameters in the top part of main function.
```sigma_sq```: migrational variance $\sigma^2$  
```s1```: selection coefficient in environment 1, $s_1$  
```s2```: selection coefficient in environment 2 (with inversed sign), $s_2$  
```r1```: the radius of the target patch in $d=1$ case  
```r2```: the radius of the target patch in $d=2$ case  
```r3```: the radius of the target patch in $d=3$ case  

2. Compile files using the following command:  
```
  g++ calculate_est.cpp -Wall -Wextra -O3 -std=c++17 -lgmp -lgmpxx -o XXX.out
```
It is required to install gmp and boost in advance.

3. Run the executive file:
```
./XXX.out
```

#### Output files
Final files: "output_int_1d.txt", "output_int_2d.txt", "output_int_3d"  
Each file contains results for corresponding spatial dimension.  

1st col: the distance between the mutation origin and the center of the target patch, $r$  
2nd col: estabilishment probability, $u_d(r)$  
3rd col: the relative contribution, $\kappa(r)$  
4th col: $\int_{x\in\text{space}} u_d(x)dx$

### "calculate_wait_time.cpp"  
#### Usage
1. This program uses outputs from "calculate_est.cpp", so please make sure that "calculate_est.cpp" was run with correct parameters before running this code.  
2. Set parameters in the top part of main function.  
```sigma_sq```: migrational variance $\sigma^2$  
```s2```: selection coefficient in environment 2 (with inversed sign), $s_2$  
```r```: radius of the target patch  
```rs```: radius of the source patch  
```rho```: population density, $\rho$  
```d_max```: maximum distance investigated between the patches  
3. Compile files using the following command:  
```
  g++ calculate_wait_time.cpp -Wall -Wextra -O3 -std=c++17 -o XXX.out
```
It is required to install boost in advance.

4. Run the executive file:
```
./XXX.out
```

#### Output file
Final file: "lambda.txt"  

1st col: distance between the source and target patch, $D$  
2nd col: $1/\lambda_{\text{mig}}$ for $d=1$  
3rd col: $T_{\text{wait}}$ for $d=1$  
4th col: $1/\lambda_{\text{mig}}$ for $d=2$  
5th col: $T_{\text{wait}}$ for $d=2$  
6th col: $1/\lambda_{\text{mig}}$ for $d=3$  
7th col: $T_{\text{wait}}$ for $d=3$  

### "calculate_wait_time_RC.cpp"
#### Usage
1. Set parameters in the top part of main function.
```sigma_sq```: migrational variance $\sigma^2$  
```s1```: selection coefficient in environment 1, $s_1$  
```s2```: selection coefficient in environment 2 (with inversed sign), $s_2$  
```r```: radius of the target patch  
```rs```: radius of the source patch  
```rho```: population density, $\rho$  
```max_d```: maximum distance investigated between the patches 


2. Compile files using the following command:  
```
  g++ calculate_wait_time_RC.cpp -Wall -Wextra -O3 -std=c++17 -o XXX.out
```

3. Run the executive file:
```
./XXX.out
```

#### Output file
Final file: "predict_RC.txt"  

1st col: distance between the source and target patch, $D$  
2nd col: $1/\lambda_{\text{mig}}$ for $d=1$  
3rd col: $1/\lambda_{\text{mig}}$ for $d=2$  
4th col: $1/\lambda_{\text{mig}}$ for $d=3$  
