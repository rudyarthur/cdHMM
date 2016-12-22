# cdHMM
Continuous and Discrete Hidden Markov Models

cdHMM is (yet another) C++ library for Hidden Markov Models (HMMs). The nice feature of cdHMM is that it allows different emission models for each hidden state. For example one can have a two state model where observations are distributed like a Gaussian in one state and like a power law in the other.

cdHMM can be used as a header only library or there is a commandline utility that wraps some of the functionality.

# Command Line

Assuming you have a relatively new version of g++ (cdHMM uses -std=c++11) typing `make` should create `hmm`.

The command line utility simply wraps some of the functionality of cdHMM and may be sufficient for your uses. e.g.
```
./hmm --in double --gaussian 1 --pareto 1 data.txt
```
will read in `data.txt` (assumed to be a single column of double precision numbers) and fit using a mixed emission model with one Gaussian state and one power law state.

Type `./hmm` without arguments for a full list of options.

Typing `make test` will make the test program, and running `./test` will give a list of options.

# Library

cdHMM can be included in your own code simply by adding
```
#include "[pathtocdHMM]/cdhmm.hpp"
...
using namespace cdHMM;
...
```
or, if you prefer to avoid global namespaces, by putting cdHMM:: in front of declarations like
```
cdHMM::gaussianHMM<double> hmm(2,0,1000);
```
cdHMM hmm objects can be constructed explicitly as above or using the HMM factory
```
HMM<double>* hmm = cdHMM::HMMFactory<double>::generateHMM("gaussian", 2, 0, 1000);
... do something with hmm ...
delete hmm;
```
The second form is useful when you want a mixed emission model. For example the two state model, where one state has Gaussian emission and the other has power law emission (called Pareto in cdHMM)
```
HMM<double>* hmm1 = cdHMM::HMMFactory<double>::generateHMM("gaussian", 1, 0, 1000);
HMM<double>* hmm2 = cdHMM::HMMFactory<double>::generateHMM("pareto", 1, 0, 1000);
vector< HMM<double>* > hmms = {hmm1, hmm2};
multiHMM<double> hmm(hmms,0,1000);
...
delete hmm1;
delete hmm2;
```
The `multiHMM` type show here allows mixed emission models.

The most useful HMM functions are

* `hmm.info()`: prints info about the particular hmm object.
* `hmm.init()`: puts hmm object in ready state to fit data
* `hmm.fit(data, eps)`: Uses EM algorithm to fit transition & emission probabilities
* `hmm.set_state_viterbi(data)`: Uses Viterbi algorithm to calculate most likely state sequence.
* `hmm.generate_seq(O, num, true)`: Generate num synthethic data points from hmm, save in O and also output to console.

And public data members 

* `hmm.A`: Transition Matrix `N x N vector<vector<double> >`
* `hmm.B`: Emission Matrix `N x M vector<vector<double> >`
* `hmm.pi`: Starting Probabilities `N vector<double>`
* `hmm.state`: State found by fit or by set_state_viterbi `num_observations vector<unsigned>`
* `hmm.maxA`: Best Transition Matrix found by fit command `N x N vector<vector<double> >`
* `hmm.maxB`: Best Emission Matrix found by fit command `N x M vector<vector<double> >`
* `hmm.maxpi`: Best Starting Probabilities found by fit command `N vector<double>`
* `hmm.maxlhood`: Maximum log likelihood found by fit command `double`

# Examples:
We can do standard discrete hmm things. For example let's take the first chapter of Pride and Prejudice (provided in
data/pridech1.txt, obtained from [project Gutenberg](https://www.gutenberg.org/ebooks/1342)), stripping out all puctuation and capitalisation:
```
./hmm --in char --tolower --discrete 2 --precision 1e-10 --restarts 10 example_data/pridech1.txt
```
```
Transition Matrix
0.260607	0.739393	
0.757925	0.242075	
Emission Matrix
i:	0.114304	0.000000
t:	0.016926	0.116475
 :	0.378303	0.029674
s:	0.000000	0.105833
a:	0.110043	0.002471
r:	0.000000	0.093968
u:	0.034273	0.023701
h:	0.000000	0.099188
n:	0.000000	0.117697
v:	0.000000	0.017560
e:	0.188348	0.000000
l:	0.006117	0.051152
y:	0.000000	0.052204
c:	0.000000	0.028950
k:	0.001918	0.009898
o:	0.129576	0.000000
w:	0.000000	0.037018
d:	0.000000	0.056950
g:	0.018788	0.014428
m:	0.000000	0.056950
p:	0.000000	0.013763
f:	0.000000	0.040814
b:	0.000000	0.024678
x:	0.000000	0.002848
j:	0.000000	0.001424
z:	0.001405	0.001406
q:	0.000000	0.000949
```
So we have managed to discover that there are two classes of character: toughly 'a', 'e', 'i', 'o' and ' ' in state 1 and everything else in state 2. Surprisingly 'u' is ambiguous - not much more likely in state 1 or state 2, but this is a relatively small sample of the English language.

We can demonstrate the use of a mixed model on some publicly available data - the [Old Faithful Geyser data](http://www.stat.cmu.edu/~larry/all-of-statistics/=data/faithful.dat). Trying to model the eruption durations we try a model with one Gaussian state and one gamma state.
```
./hmm --in double --col 1 --gaussian 1 --gamma 1 --precision 1e-10 --restarts 10 --generate 1000 example_data/oldfaithful.txt
```
![alt text](https://github.com/rudyarthur/cdHMM/blob/master/example_data/oldfaithful.png)


