<h1 align='center'>Code for the paper "A binomial-adjusted θ-trapezoidal tau-leaping method 
for chemical kinetic systems"</h1>
<h2 align='center'>by Maximilian Scott, Sam McCallum and James Foster</h2>


The paper is available at [insert link].

This repository contains code for simulating chemical reaction systems using tau-leaping methodologies, in particular the BθTTLe numerical scheme.
For example consider the standard tau-leaping, θ-trapezoidal and BθTTLe methods applied to a kinetic 
model with 11 reacting species and 11 reaction channels. We can compare their respective endpoints to 
those of an exact sampling methodology using Fréchet Inception Distance vs the bulk of the associated 
computational cost.

<p align="middle">
    <img src="./experiments/experiment_1_FID.png" width=45%>
</p>

Additional experimental data can be found `data`. 


