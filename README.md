## Homogenisation for arbitrary arrangements of obstructions in cardiac fibrosis

This MATLAB code accompanies the journal article "Homogenisation for the monodomain model in the presence of microscopic fibroticstructures",
by Brodie Lawson, Rodrigo Weber dos Santos, Ian Turner, Alfonso Bueno-Orovio, Pamela Burrage and Kevin Burrage. All code created by Brodie Lawson,
in co-operation with the other authors.

This code includes:
  1. A numerical solver for the monodomain model on a regular 2D grid, allowing for heterogeneities in conduction and/or cell model
  2. Methods for calculation of effective conductivity tensors through homogenisation by volume averaging
  3. A visualisation tool
  

Numerical solution of the monodomain model is achieved by using a control volume finite element approach, and the second-order generalisation of
Rush-Larsen timestepping published by Perego and Veneziani (2009). Homogenisation by volume averaging is credit to Whitaker (1999) and the
homogenisation approach allows for non-periodic boundary conditions discussed for example by Wu et al. (2002). Full reference information is provided
in the accompanying journal article.
