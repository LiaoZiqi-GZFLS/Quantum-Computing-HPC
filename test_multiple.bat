@echo off

setlocal enabledelayedexpansion

:input
set /p "filePath=Please input the path of cpp: "

if exist "!filePath!" (
    echo File exists!
) else (
    echo File doesn't exist!
    goto :input
)

set  "remotePath=/work/sustcsc_02/workspace"
set  "remotePath2=/work/sustcsc_02/workspace/code"
set  "user=sustcsc_02"
set  "password=yeshenhaike"
set  "port=18188"
set  "sever=172.18.6.40"

set "localFilePath=file"

echo %path% | findstr /i "putty"  >nul
if %errorlevel% equ 0 (
    echo PuTTY exists
) else (
    echo PuTTY doesn't exist
	:: 打开 PuTTY 下载网站
    	echo try at https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html
	start https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html
	pause>nul
	exit
)

echo Submitting...
pscp -ssh -pw %password% -P %port% "!filePath!" %user%@%sever%:"!remotePath2!"
if %errorlevel% equ 0 (
        echo Succeed!
    ) else (
        echo Fail.
	pause>nul
	exit
    )

set /p "input=Test 22.3GB demo in computing machine or doing simple test?(y/n)  "
if /i "%input%" == "y" (
	plink -ssh -pw %password% -P %port% %user%@%sever%  -batch "%remotePath%/auto_submit.sh auto_run.sh"
	if %errorlevel% equ 0 (
        echo Succeed!
   	 ) else (
        echo Fail.
		pause>nul
		exit
    )

	echo Default:
	plink -ssh -pw %password% -P %port% %user%@%sever%  -batch "cat %remotePath%/ref_domo_result.txt"
) else (
    echo Testing...
    plink -ssh -pw %password% -P %port% %user%@%sever%  -batch "%remotePath%/auto_submit.sh auto_test.sh"
    if %errorlevel% equ 0 (
            echo Succeed!
    ) else (
            echo Fail.
        pause>nul
        exit
    )
)

endlocal

pause>nul