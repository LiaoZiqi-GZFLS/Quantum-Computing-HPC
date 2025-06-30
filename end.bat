@echo off
plink -ssh -pw yeshenhaike -P 18188 sustcsc_02@172.18.6.40  -batch "/work/sustcsc_02/workspace/stop.sh"
pause>nul
exit /b %errorlevel%