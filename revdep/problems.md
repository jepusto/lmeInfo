# scdhlm

<details>

* Version: 0.5.2
* GitHub: https://github.com/jepusto/scdhlm
* Source code: https://github.com/cran/scdhlm
* Date/Publication: 2021-01-07 19:30:02 UTC
* Number of recursive dependencies: 93

Run `revdep_details(, "scdhlm")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in 'scdhlm-Ex.R' failed
    The error most likely occurred in:
    
    > ### Name: CI_g
    > ### Title: Calculates a confidence interval for a standardized mean
    > ###   difference effect size
    > ### Aliases: CI_g
    > 
    > ### ** Examples
    > 
    ...
    [1] 0.9143684 2.0046719
    > 
    > Laski_HPS <- with(Laski, effect_size_MB(outcome, treatment, case, time))
    > CI_g(Laski_HPS, symmetric = FALSE)
    [1] 0.8786207 2.0639828
    > 
    > Laski_g_mlm <- g_mlm(Laski_RML, p_const = c(0,1), r_const = c(1,0,1), returnModel = TRUE)
    Error in g_mlm(Laski_RML, p_const = c(0, 1), r_const = c(1, 0, 1), returnModel = TRUE) : 
      unused argument (returnModel = TRUE)
    Execution halted
    ```

*   checking tests ...
    ```
      Running 'testthat.R'
     ERROR
    Running the tests in 'tests/testthat.R' failed.
    Last 13 lines of output:
      > 
      > test_check("scdhlm")
      [ FAIL 1 | WARN 0 | SKIP 8 | PASS 89 ]
      
      == Skipped tests ===============================================================
      * Auxiliary dataset not included in package. (1)
      * empty test (7)
      
      == Failed tests ================================================================
      -- Error (test-g_mlm.R:22:3): g_mlm() is imported appropriately. ---------------
      Error in `g_mlm(Laski_RML1, p_const = c(0, 1), r_const = c(1, 0, 1), returnModel = TRUE)`: unused argument (returnModel = TRUE)
      
      [ FAIL 1 | WARN 0 | SKIP 8 | PASS 89 ]
      Error: Test failures
      Execution halted
    ```

## In both

*   checking whether package 'scdhlm' can be installed ... WARNING
    ```
    Found the following significant warnings:
      Warning: package 'nlme' was built under R version 4.1.3
    See 'D:/Amanda/Pusto/git/lmeInfo/revdep/checks/scdhlm/new/scdhlm.Rcheck/00install.out' for details.
    ```

