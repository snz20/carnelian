******************
Running unit tests
******************
Command:
  python basictest_carnelian.py

Desired effect:
  If run successfully, you will see the following message at the end of the console output:
  ----------------------------------------------------------------------
  Ran 32 tests in XX.XXXs

  OK

Note:
- The above script will attempt to run the tests using 5 cpus. If you do not have access to multiple processors, run the 
  following command:
    python basictest_carnelian.py -n 1
- The above script assumes that you have FragGeneScan installed in the util/ext/FragGeneScan directory under the carnelian's
  source directory and you have a working vowpal-wabbit and R installation.

**********************
Running advanced tests
**********************
Command:
  python advancedtest_carnelian.py

Desired effect:
  If run successfully, you will see the following message at the end of the console output:
  ----------------------------------------------------------------------
  Ran 2 tests in XX.XXXs

  OK

Note:
- The above script will attempt to run the tests using 5 cpus. If you do not have access to multiple processors, run the 
  following command:
    python advancedtest_carnelian.py -n 1
- The above script assumes that you have FragGeneScan installed in the util/ext/FragGeneScan directory under the carnelian's
  source directory and you have a working vowpal-wabbit and R installation.
