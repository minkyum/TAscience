#!/bin/bash
echo Submitting tile $1
echo Processing year $2
R --vanilla < /usr3/graduate/mkmoon/GitHub/TAscience/mcd12q2/MCD12Q2C6_AnnualPhenology.R --args -tile $1 -year $2

