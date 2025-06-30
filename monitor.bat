@echo off
:_book
plink -ssh -pw yeshenhaike -P 18188 sustcsc_02@172.18.6.40  -batch "mpstat -P ALL 1 1"
goto _book