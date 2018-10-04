

#!/bin/bash
iter=1
iterMax=10
myExe=dmft

while [  $iter -lt $iterMax ];
    do
       $myExe params $iter
       iter=$[iter+1]
    done


