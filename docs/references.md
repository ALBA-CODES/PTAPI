# References

What each paper is used for in the API. The bibliographic data of the entries
marked with † were written from memory and should be checked against the
publisher before being cited.

## Parallel Tempering and the temperature ladder

- GEYER, C. J. Markov chain Monte Carlo maximum likelihood. In: *Computing Science and Statistics: Proceedings of the 23rd Symposium on the Interface*, 1991. — The method.
- HUKUSHIMA, K.; NEMOTO, K. Exchange Monte Carlo method and application to spin glass simulations. *Journal of the Physical Society of Japan*, v. 65, n. 6, p. 1604–1608, 1996. — The method.
- HUKUSHIMA, K. Domain-wall free energy of spin-glass models: numerical method and boundary conditions. *Physical Review E*, v. 60, p. 3606–3613, 1999. — Update type 2 (equal exchange rates).
- KATZGRABER, H. G.; TREBST, S.; HUSE, D. A.; TROYER, M. Feedback-optimized parallel tempering Monte Carlo. *Journal of Statistical Mechanics: Theory and Experiment*, P03018, 2006. — Update types 3 and 4, round trips.
- ROBERTS, G. O.; GELMAN, A.; GILKS, W. R. Weak convergence and optimal scaling of random walk Metropolis algorithms. *The Annals of Applied Probability*, v. 7, n. 1, p. 110–120, 1997. — The 23.4% rate.
- ATCHADÉ, Y. F.; ROBERTS, G. O.; ROSENTHAL, J. S. Towards optimal scaling of Metropolis-coupled Markov chain Monte Carlo. *Statistics and Computing*, v. 21, p. 555–568, 2011. — Exchange rate of 23% (update type 1, default of `calibrateLadder`).
- MIASOJEDOW, B.; MOULINES, E.; VIHOLA, M. An adaptive parallel tempering algorithm. *Journal of Computational and Graphical Statistics*, v. 22, n. 3, p. 649–664, 2013. arXiv:1205.1076. † — Update type 5 and the correction rule of `calibrateLadder`.
- SYED, S.; BOUCHARD-CÔTÉ, A.; DELIGIANNIDIS, G.; DOUCET, A. Non-reversible parallel tempering: a scalable highly parallel MCMC scheme. *Journal of the Royal Statistical Society: Series B*, v. 84, n. 2, p. 321–350, 2022. arXiv:1905.02939. † — Update type 6, even-odd exchange scheme, communication barrier, doubling of the update interval.
- OKABE, T.; KAWATA, M.; OKAMOTO, Y.; MIKAMI, M. Replica-exchange Monte Carlo method for the isobaric–isothermal ensemble. *Chemical Physics Letters*, v. 335, p. 435–439, 2001. † — Even-odd exchange scheme.
- BEN-AMEUR, W. Computing the initial temperature of simulated annealing. *Computational Optimization and Applications*, v. 29, n. 3, p. 369–385, 2004. † — Idea related to `calibrateTemperatures`.

## Experiments and statistics

- MANN, H. B.; WHITNEY, D. R. On a test of whether one of two random variables is stochastically larger than the other. *The Annals of Mathematical Statistics*, v. 18, n. 1, p. 50–60, 1947. †
- WILCOXON, F. Individual comparisons by ranking methods. *Biometrics Bulletin*, v. 1, n. 6, p. 80–83, 1945. †
- FRIEDMAN, M. The use of ranks to avoid the assumption of normality implicit in the analysis of variance. *Journal of the American Statistical Association*, v. 32, n. 200, p. 675–701, 1937. †
- HOLM, S. A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics*, v. 6, n. 2, p. 65–70, 1979. †
- DEMŠAR, J. Statistical comparisons of classifiers over multiple data sets. *Journal of Machine Learning Research*, v. 7, p. 1–30, 2006. † — The methodology of `PTComparison`.
- AIEX, R. M.; RESENDE, M. G. C.; RIBEIRO, C. C. TTT plots: a perl program to create time-to-target plots. *Optimization Letters*, v. 1, p. 355–366, 2007. † — `time_to_target.csv`.
- LÓPEZ-IBÁÑEZ, M.; DUBOIS-LACOSTE, J.; PÉREZ CÁCERES, L.; BIRATTARI, M.; STÜTZLE, T. The irace package: iterated racing for automatic algorithm configuration. *Operations Research Perspectives*, v. 3, p. 43–58, 2016. † — irace.

## Read, not implemented

- ŁĄCKI, M. K.; MIASOJEDOW, B. State-dependent swap strategies and automatic reduction of number of temperatures in adaptive parallel tempering algorithm. *Statistics and Computing*, v. 26, p. 951–964, 2016. arXiv:1406.3218. †
- SYED, S.; ROMANIELLO, V.; CAMPBELL, T.; BOUCHARD-CÔTÉ, A. Parallel tempering on optimized paths. *ICML*, 2021. arXiv:2102.07720. †
- DELACOUR, C.; SAJEEB, M. M. H.; HESPANHA, J. P.; CAMSARI, K. Y. Two-dimensional parallel tempering for constrained optimization. arXiv:2506.14781, 2025.
- SIM, K.; CHOI, S.; PARK, J. Solving integer linear programming with parallel tempering. arXiv:2605.29366, 2026.
- MIYATA, T.; ARAI, S.; TAKABE, S. Refined gradient-based temperature optimization for the replica-exchange Monte-Carlo method. arXiv:2601.13542, 2026.

## The API

- ALMEIDA, A. L. B. *Revisitando o revenimento paralelo: computação de alto desempenho e aplicação em pesquisa operacional*. Doctoral thesis, Universidade Federal de Ouro Preto, 2024. <https://www.repositorio.ufop.br/handle/123456789/18678> — The dataflow implementation for multi-core CPUs (chapter 5) and the three temperature updates of version 1. See [`CITATION.cff`](../CITATION.cff).
