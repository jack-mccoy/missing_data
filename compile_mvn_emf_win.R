# Compile mvn_emf.f90 for Windows
# Assumes Rtools 4.2 or more recent is installed

#   this is poorly documented, but see here: 
#   https://stackoverflow.com/questions/24044169/build-shared-object-from-fortran-source-with-packagee-g-lapack-for-r
#   and here:
#   https://cran.r-project.org/bin/windows/Rtools/rtools40.html

# create make file
writeLines(
  'PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)', con = 'makevars'
)

# clean
system('rm mvn_emf_win.so')

# run compile command
system('R CMD SHLIB mvn_emf.f90 -o mvn_emf_win.so')

# clean up
file.remove('makevars')
