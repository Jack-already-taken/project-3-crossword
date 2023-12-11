Link to the github repository: https://github.com/eec289q-f23/project-3-crossword-compiler-Jack-already-taken.git

## Solve Order Optimization
I started off creating a function that can generate a more optimized solve order so that the search program can find the solution faster without needing to traverse though all possible ways to fill the grid. I started by creating a dynamic fill order function ```findNextFill()``` that attempts to find the entry on the grid that is crossed by most number of filled words. The filled words are composed of the seeds and the possible solutions that are filled after each recursive call of ```solve()```. When I realized that this dynamic search order is generating deterministic search entry at each recursion depth of ```solve()```, I created a static search order function ```filledAwareOrder()``` that repeatedly calls ```findNextFill()``` and update the choice made by the dynamic search function on a simulated grid until empty entries on the grid are exhausted.

The code with optimized search order function have the best performance speedup when the seeds are placed near the center of the grid. For example:
```console
$ srun -w agate-1 --reservation=eec289q -t 5 --cpu-bind=ldoms -c 1 ./xword -m 50 -g grids/nyt-230916.txt -w wordlists/spreadthewordlist.txt -s 0,1,d,onahigh -s 0,13,d,lilbaby -s 0,14,d,sclass -s 1,8,d,goingtoetotoe -s 6,1,a,herdmentality -s 8,1,a,whatsthepoint
```
The above test case took the code with optimized search order 0.993868 seconds while it took the original code 5.19568 seconds to run. However, for the benchmarking case of the leaderboard:
```console
$ srun -w agate-1 --reservation=eec289q -t 3 --cpu-bind=ldoms -c 16 ./xword -g grids/nyt-230916.txt -w wordlists/spreadthewordlist.txt -s 0,0,a,scales -m 40
```
Since the seed is at the top left corner, the optimized search order code doesn't have any performance increase over the original code.

## Parallelizing solve()
#### Issue with return statements
After chaning the for loop in the function to cilk_for, I dealt with this issue by using an atomic boolean flag which is declared as the member of the ```Solve``` class. The flag is initialized to false and is only set to true through compare-and-swap in the condition ```i == solveOrder.size()```, whihc means when the solve has found the solution. I added a condition to check the value of the flag at the beginning of the cilk_for loop, and if the flag is true, the rest of the cilk_for iteration will be continued for a faster exit.

#### Race condition with Wordlist solution
Since the cilk_for loop involves assigning the ```solution``` wordlist at the same search position with differnt words from the searchspace, this becomes a race condition that could corrupt the solution from valid to invalid. Therefore, I modified the code structure of ```solve()``` such that it takes an extra argument, which is the copy of the ```solution``` Wordlist each time when it is recursively called. When reaching the part of code that changes the flag from false to true, the ```solve()``` function will copy the content of the ```solution``` Wordlist into the class member variable ```finalSolution```. Since the compare-and-swap function is guaranteed to assign the flag only once during the entire code execution, the ```finalSolution``` Wordlist will be free from datarace. However, the solution that the ```solve()``` function generates will no longer be deterministic, since it will end up at different correct solutions based on how the scheduler schedules the threads.

I found that adding only cilk_for parallelization, which means full breadth-first search, results in the fastest code. This is because I tried changing the code such that at certain recursion depth of ```solve()```, it will do depth-first search instead of breadth-first search, meaning that it will use regular for loop when the recursion depth is deep enough. This resulted in a significant slowdown.

## Parallelizing computeNewSearchSpaceGivenGuess()
I parallelized this function mainly because when I checked the perf report, I found that the function ```constrainWordlist()``` is the bottleneck of the program. Since the function ```computeNewSearchSpaceGivenGuess()``` calls this function in a for loop, I parallelized the function ```computeNewSearchSpaceGivenGuess()``` using cilk_for. This doesn't result in significant speedup over the best case runtime, but it makes the runtime of the program more consistently stay at the best case with different inputs.

## Resulting Speedup
For the leaderboard benchmarking case, I was able to get this runtime:
```
Elapsed Time: 0.756549 seconds
Solution found.
Grid is (15, 15) and VALID
scales###brians
taliban#saidboo
rcadome#hirsute
ehs#narrate#bsa
eeka#rvers#gaps
ttada#oop#aokay
#snoreupastorm#
###ramsesiii###
#chickenandegg#
abone#nik#esale
grog#benny#trad
las#borgias#ari
editing#fridges
tierney#entreat
sorest###saidto
```
For the test case:
```console
$ srun -w agate-1 --reservation=eec289q -t 1 --cpu-bind=ldoms -c 16 ./xword -m 50 -g grids/nyt-230916.txt -w wordlists/spreadthewordlist.txt -s 5,0,a,aggie
```
which the original program took hours to finish, the optimized code can obtain this runtime:
```
Elapsed Time: 15.422 seconds
Solution found.
Grid is (15, 15) and VALID
bating###allies
tvidols#lieidle
redspot#imadeit
ari#evicted#act
tank#elote#amie
aggie#ell#jdate
#espritdecorps#
###pantsrole###
#chessopenings#
blare#hed#eared
lots#melba#leno
use#dyelots#ade
russell#orestes
troikas#killers
sensei###alonso
```

## Potential Issue
The issue I have with the current code is that cilksan will always say that my current code has race condition even though my program can reliably find the correct solution for all the test cases. Example below:
```console
...
Elapsed Time: 5.07645 seconds
Solution found.
Grid is (15, 15) and VALID
cosmos###broils
anoints#graphic
baldcap#outsell
aha#egalite#aba
nina#stone#iras
aggie#rag#cribs
#herdmentality#
###didasolid###
#whatsthepoint#
broth#mat#subin
aime#perot#many
rte#tanktop#sim
restart#onestep
emigres#eyelash
denied###skorts


Cilksan detected 22 distinct races.
Cilksan suppressed 88799 duplicate race reports.
```
Because of this, I also tried implementing ```solution``` reducer in the file ```main_reducer.cpp```. The code executes but faces a race condition that affects the ability of the program to find the correct solution. If time permits I would really like to understand how my reducer implementation was not working, but please use ```main.cpp``` for grading.