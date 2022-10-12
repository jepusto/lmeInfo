# merDeriv

<details>

* Version: 0.2-4
* GitHub: https://github.com/nctingwang/merDeriv
* Source code: https://github.com/cran/merDeriv
* Date/Publication: 2022-03-11 23:30:05 UTC
* Number of recursive dependencies: 64

Run `revdep_details(, "merDeriv")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running 'tinytest.R'
     ERROR
    Running the tests in 'tests/tinytest.R' failed.
    Last 13 lines of output:
      test_merDeriv.R...............   19 tests [0;32mOK[0m 
      test_merDeriv.R...............   19 tests [0;32mOK[0m 
      test_merDeriv.R...............   19 tests [0;32mOK[0m 
      test_merDeriv.R...............   20 tests [0;31m1 fails[0m 
      test_merDeriv.R...............   20 tests [0;31m1 fails[0m 
      test_merDeriv.R...............   20 tests [0;31m1 fails[0m 
      test_merDeriv.R...............   20 tests [0;31m1 fails[0m 
      test_merDeriv.R...............   21 tests [0;31m1 fails[0m 
      test_merDeriv.R...............   22 tests [0;31m1 fails[0m 
      test_merDeriv.R...............   23 tests [0;31m1 fails[0m [0;34m1.4s[0m
      ----- FAILED[data]: test_merDeriv.R<153--153>
       call| expect_true(all(abs(nlmeres - lme4res) < 1))
       diff| Expected TRUE, got FALSE
      Error: 1 out of 23 tests failed
      Execution halted
    ```

