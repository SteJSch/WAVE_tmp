# Standard
Ucomp = "gfortran -std=legacy -c -g -cpp -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -fno-automatic -ffixed-line-length-none -finit-local-zero -funroll-loops "

#Standard, but suppress warnings
Ucomp_nowarn = "gfortran -w -std=legacy -c -g -cpp -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -fno-automatic -ffixed-line-length-none -finit-local-zero -funroll-loops "

#Standard, but with all checks at run-time
Ucomp_all = "gfortran -std=legacy -c -g -cpp -fcheck=all -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -fno-automatic -ffixed-line-length-none -finit-local-zero -funroll-loops "

#Standard, but with all checks at run-time and OpenMP (i.e. parallel processing)
Ucomp_omp = "gfortran -std=legacy -c -g -cpp -finit-local-zero -fcheck=all -fopenmp -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -ffixed-line-length-none -funroll-loops "
