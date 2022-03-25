
$onecho > cost.txt
par=i_cost                       rng=cost!b4:h13
$offecho

$CALL GDXXRW.exe    "in.xlsx"           @cost.txt

*$GDXIN "in.xlsx"
*$LOADdc i_cost

* ---



$Onecho > task.txt
par=i_load     rng=Sheet1!A1:D11    rdim=1  cdim=1
$Offecho

$call GDXXRW.exe load.xlsx  @task.txt

* ---

$Onecho > task.txt
par=avail     rng=Sheet1!A1:G13    rdim=1  cdim=2
$Offecho

$call GDXXRW.exe avail.xlsx  @task.txt


 