#!/bin/bash


echo Submitting $1
R --vanilla < /usr3/graduate/mkmoon/GitHub/MCD12Q2/get_MCD43A2.R $1

echo Submitting $1
R --vanilla < /usr3/graduate/mkmoon/GitHub/MCD12Q2/get_MCD43A4.R $1


