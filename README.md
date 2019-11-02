Author
-------

Hanjie Pan
Audiovisual Communications Laboratory ([LCAV](http://lcav.epfl.ch)) at [EPFL](http://www.epfl.ch).<br>

<img src="./html/LCAV_anim_200.gif">

#### Contact

[Hanjie Pan](mailto:panhanjie[at]gmail[dot]com)<br>
EPFL-IC-LCAV<br>
BC 322 (Bâtiment BC)<br>
Station 14<br>
1015 Lausanne

Abstract
--------

Image up-sampling aims at reconstructing a high resolution image from a single low resolution one. 
It is essential to exploit prior knowledge on the reconstructed image in order to better condition 
this severely ill-conditioned inverse problem. We present a novel edge modelling framework based on 
a spatial domain interpretation of annihilation of curves with finite rate of innovation. 
More specifically, we define a continuous domain mask function that vanishes around large gradients, 
i.e., around edges. The mask function is reconstructed by minimising its product with the gradient image. 
We show that accurate edge models are reconstructed by assuming a simple local linear edge model. 
Based on the same idea, we further combine these local edge models and build a global one, 
which serves as an edge-preserving constraint in image up-samplings. Moreover, we propose 
an efficient alternating direction method of multipliers (ADMM) to solve the up-sampling 
problem numerically. Experiments with natural images demonstrate the effectiveness 
of the global edge model in improving the quality of the up-sampled images thus achieving state 
of the art performance.

Package Description
-------------------

This repository contains the code to reproduce the results of Chapter 6 of [*Looking beyond Pixels Theory, Algorithms and Applications of Continuous Sparse Recovery*](https://infoscience.epfl.ch/record/255725). It contains a Matlab implementation of the proposed algorithm.

Run image up-sampling
--------------------------------------
    % image up-sampling with known ground-truth
    main_upSamp.m
    
    % blind image up-sampling
    main_blindUpSamp.m

License
-------

Copyright (c) 2018, Hanjie Pan<br>
The source code is released under the [MIT](LICENSE.txt) license.