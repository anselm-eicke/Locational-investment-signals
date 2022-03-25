
*-----------------------------------
* LOAD INPUT PARAMETERS FROM EXCEL

$onecho > cost.txt
par=i_cost                       rng=cost!b4:h13
$offecho

* Convert XLSX to GDX file and load
$if not set LOADDATA   $set LOADDATA  1

$CALL GDXXRW.exe    "in.xlsx"           @cost.txt
