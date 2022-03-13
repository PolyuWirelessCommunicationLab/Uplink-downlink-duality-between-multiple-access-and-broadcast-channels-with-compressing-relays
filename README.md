# Uplink-downlink-duality-between-multiple-access-and-broadcast-channels-with-compressing-relays
This code is for paper: [L. Liu, Y-F. Liu, P. Patil, and W. Yu, "Uplink-downlink duality between multiple-access and broadcast channels with compressing relays," IEEE Trans. Inf. Theory, vol. 67, no. 11, pp. 7304-7337, Nov. 2021.](https://ieeexplore.ieee.org/abstract/document/9514556/)
# Abstract
Uplink-downlink duality refers to the fact that under a sum-power constraint, the capacity regions of a Gaussian multiple-access channel and a Gaussian broadcast channel with Hermitian transposed channel matrices are identical. This paper generalizes this result to a cooperative cellular network, in which remote access-points are deployed as relays in serving the users under the coordination of a central processor (CP). In this model, the users and the relays are connected over noisy wireless links, while the relays and the CP are connected over noiseless but rate-limited fronthaul links. Based on a Lagrangian technique, this paper establishes a duality relationship between such a multiple-access relay channel and broadcast relay channel, under the assumption that the relays use compression-based strategies. Specifically, we show that under the same total transmit power constraint and individual fronthaul rate constraints, the achievable rate regions of the Gaussian multiple-access and broadcast relay channels are identical, when either independent compression or Wyner-Ziv and multivariate compression strategies are used. The key observations are that if the beamforming vectors at the relays are fixed, the sum-power minimization problems under the achievable rate and fronthaul constraints in both the uplink and the downlink can be transformed into either a linear programming or a semidefinite programming problem depending on the compression technique, and that the uplink and downlink problems are Lagrangian duals of each other. Moreover, the dual variables corresponding to the downlink rate constraints become the uplink powers; the dual variables corresponding to the downlink fronthaul constraints become the uplink quantization noises. This duality relationship enables an efficient algorithm for optimizing the downlink transmission and relaying strategies based on the uplink.
# Readme
ul_vs_dl_wcmc_versus_ind.m is the main function.
# Citation
@article{liu2021duality,<br> 
  title={Uplink-downlink duality between multiple-access and broadcast channels with compressing relays},<br> 
  author={L. Liu, Y-F. Liu, P. Patil, and W. Yu},<br> 
  journal={IEEE Trans. Inf. Theory},<br> 
  volume={67},<br> 
  number={11},<br> 
  pages={7304--7337},<br> 
  month={Nov.},<br>
  year={2021},<br> 
  publisher={IEEE}<br> 
}<br> 
# Note
The code is provided for the benefit of better understanding the results, and is not meant to be used in production.
