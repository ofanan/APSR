# APSR
This projects implements a part of the APSR (Adaptive Partial State Random) algorithm for VM placement with high throughput and provalbe success guarantees. The APSR algorithm is described in the paper:

I. Cohen, G. Einziger, M. Goldstein, Y. Sa’ar, G. Scalosub, and E. Waisbard. [Parallel VM Deployment with Provable Guarantees](https://www.researchgate.net/profile/Itamar-Cohen-2/publication/351449290_Parallel_VM_Deployment_with_Provable_Guarantees/links/6098afaaa6fdccaebd1d82f5/Parallel-VM-Deployment-with-Provable-Guarantees.pdf), IFIP Networking, 2021, pp. 1-9.

The project includes an implementation of the _MaximizeParallelism_ (Algorithm 2 in the above paper), that maximize the number of parallel VM placement schedulers, while providing success guarantees. The algorithm is implemented both in Java, and in Python.
