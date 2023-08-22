---
layout: default
title: Coding conventions and regression testing
tagline: GreeX Testing
description: ' '
---
# Coding conventions

 We utilize modern Fortran features such as object-oriented programming and intrinsic procedures that are available in Fortran 2008. We have developed a clear interface between the module code (our library) and the client code (many-body perturbation theory code), promoting better modularity and reusability. Additionally, we use allocatable arrays and automatic deallocation to simplify the code and avoid memory-related issues, such as leaks and dangling pointers. The implementation is robust and reliable, as we use error handling techniques to highlight and recover from exceptions. 
 
# Regression testing

To ensure the reproducibility of our results, we have implemented regression tests using the pytest infrastructure, combined with CMake and CTest. The regression tests rigorously verify the accuracy of the frequency grids and integration weights and the forward-backward cosine error. Further details on the library structure, contents, build system, and testing are available in the various README.md files at the top level of the repository and the individual component directories. 

<button onclick="goBack()">Go Back</button>

<script>
function goBack() {
  window.history.back();
}
</script>
