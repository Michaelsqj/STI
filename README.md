# STI
**Acknowledgments to Prof. Xu Li from Radiology JHU**

**Modified based on   `MMSR-STI`**

## Use guide
Run code in `main.m`

Three parts: `phantom`, `forward`, `inverse`

## Development Plan
- Mar. 10 

    Build the framework

    write, run, visualize phantom

    write, run, visualize forward

    modified `conv_kernel_tensor` -> `conv_kernel_tensor_arr`

- Mar. 13

    Finished revising `STI_inverse`

    **TODO**

    run `phantom`

- Mar. 14

    Finished revising `STI_forward` and `gen_phantom`

    **TODO**

    run `phantom`
    
- Mar. 19
    
    Finished 'STI_forward'
    
    **TODO**
    
    Debug `STI_inverse` in `MMSR-STI` and `STI`

- Mar. 20

    Tested 'STI_inverse', temporarily reduce the wG regularization.

    Finished 'evaluate'

    **Problem**

    - Eig slows down on reconstructed chi.

    - the percentage of wG is not easy to determine. in simulation, almost 98% of
chiavg(QSM) is 0; data from real world might be different.

- Mar. 22

    Run `MMSR-STI`, added chemical shift

    **TODO**

    Run the inverse for enough iters (e.g., 1000), and compare the results, between TV
regularization, 9 elements tensor, and chemical shift.
