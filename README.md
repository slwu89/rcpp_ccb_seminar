# CCB Skills Seminar: Combining C++ and R with Rcpp

May 9, 2019: Sean Wu (slwu89@berkeley.edu)

## Directory Structure

### standalone

This directory uses the running example of a stochastic SIR model to demonstrate several concepts in linking R and C++.
 
  * foreach: contains a small example of pitfalls when running Rcpp-exported C++ code in parallel and how to fix it
  * gillespie-SIR-R: contains the full-R example of the simulation model
  * gilespie-SIR-CXX: contains the single-file C++ version of the simulation model
  * gillespie-SIR-CXX-headers: contains the C++ version of the simulation model where a RNG class is stored in a seperate header file and linked to the main file
  
### gillespieCXX

This is an example R package which uses Rcpp to integrate C++ code. It illustrates several concepts including breaking the code into header (.h) and implementation (.cpp), as well as returning STL objects to R with `Rcpp::wrap`, and using smart pointers in C++ code.
