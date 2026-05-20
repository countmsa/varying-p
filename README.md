The code for the results in the paper 'Assessing count measurement systems with varying probabilities of detection'. The full simulation study takes a long time to run, thus we have provided an example of one simulation run in the folder 'simRun1'. The code for the full study can be found in the folder 'sim study'. The code to run the data example can be found in the folder 'Bubble5Fit'. 


##Simulation Run
We provide an example simulation run as this will run quickly. The optimization for all plans will run in under 30 seconds and the Fisher informations will run in approximately 30 seconds - 1 minute. In the simulation study, we perform 1000 replicates for each combination. To run one replication of combination $(\lambda, \theta, \alpha, \beta)$ = $(3,0.1,9,3)$ for $s=50$ and $r=4$, access the file "sim1run.R" in the folder "sim1run". The data and necessary source files are also in the folder "sim1run". After running this file, the results can be found in the following objects:

'mle' a 5x4 matrix where each row is a plan in the order of Y, NF, F, N, C, and the columns are the parameters in order of $\lambda, \theta, \alpha, \beta$.

'asym.sd' a 5x4 matrix of the same structure as 'mles' but the elements are the asymptotic standard deviations of the corresponding parameter.

'tau.sd' a 5 element vector where each element is the asymptotic standard deviation of $\tau$ for the corresponding plan.

'loglik' a 5x2 matrix of the log likelihood values. The first column corresponds to plans Y, NF, F, N, and C. The second column is the corresponding analogous plan in the fixed p case, where Y is left as NA.

'time'  is a 5 element vector for the computational time of each plan.

##Simulation Study
The code to run the simulation study can be found in the folder "sim study". To run every combination would take over a week. The necessary function are provided, however, we do not provide the necessary precomputed sets, but do provide the functions to create them. 


##Bubble Soap Counting Example
The R scripts for the data example are in the folder "bubble5" as "Bubble5run.R" and "Bubble5runFisher.R". Due to their size (1.23GB and 1.05 GB), the necessary precomputed files can be generated or provided upon request.
The optimization takes approximately 15-25 minutes to run. The Fisher information wil take a few hours to run. The fitted results are saved in "bubble5fit.Rdata" for comparison.
