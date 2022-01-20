
*-----------------------------------
* LOAD INPUT PARAMETERS FROM EXCEL

$onecho > input.txt
par=i_cost                       rng=cost!b4:h13
par=i_load                       rng=timeseries!b3:d19                 rdim=1 cdim=1
par=i_avail                      rng=timeseries!e2:i19                 rdim=1 cdim=2
$offecho

* Convert XLSX to GDX file and load
$if not set LOADDATA   $set LOADDATA  1

$CALL GDXXRW.exe    "in.xlsx"        @input.txt