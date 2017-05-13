# Coding note
## 1d h-Based Richards equation
$$
C(h) \frac{\partial h }{\partial t}-
      \frac{\partial }{\partial z}[k\frac{\partial h }{\partial z}]
      +\frac{\partial k }{\partial z}=0
$$
where
$C$, $k$ are nonlinear function of $h$



### finite difference approximation
Using a 1st order Taylor expansion



### a Picard iteration FDM scheme
$$
C_i^{n+1,m} \times  \frac{ h_i^{n+1,m}- h_i^n }{\Delta t}-
[\frac{ k^{n+1,m}_{i+1/2}(h_{i+1}^{n+1,m+1}-h_i^{n+1,m+1})
-k^{n+1,m}_{i-1/2}(h_{i}^{n+1,m+1}-h_{i-1}^{n+1,m+1})}
{ \Delta z^2 }]
\\
+\frac{k^{n+1,m}_{i+1/2}-k^{n+1,m}_{i-1/2}}{\Delta z}
$$

### boundary condition

### for inner nodes i (where neighbours are free nodes)
$$
w_{i-1}h_{i-1}+w_{i}h_{i}+w_{i-1}h_{i-1}=b_i
$$
where the weight could be calculated from above equation easily.

### Assemble
we now have linear equation to be solve at each iteration m
$Ah=b$. A is a sparse band matrix.
Let this describes all node including those the values are known.

Define a picking up matrix $P$ that would pick up the row corresponding to free node.

we have $\hat{P}+P=I$

$$
(\hat{P}+P)Ah=b
$$

$$
PAh=b-\hat{P}Ah
$$
Sine the zero column of $P/\hat{P}$
