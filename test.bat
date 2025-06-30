@echo off

setlocal enabledelayedexpansion

:input
set /p "filePath=Please input the path of cpp: "
set "newFile=simulate.cpp"

if exist "!filePath!" (
    echo File exists!
) else (
    echo File doesn't exist!
    goto :input
)

if exist "!newFile!" (
	del "!newFile!" /F /Q
)

copy  "!filePath!" "!newFile!" /Y

set  "remotePath=/work/sustcsc_02/workspace"
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
pscp -ssh -pw %password% -P %port% "!newFile!" %user%@%sever%:"!remotePath!"
if %errorlevel% equ 0 (
        echo Succeed!
    ) else (
        echo Fail.
	pause>nul
	exit
    )

echo Connecting...
plink -ssh -pw %password% -P %port% %user%@%sever%  -batch "%remotePath%/do_test.sh"
if %errorlevel% equ 0 (
        echo Succeed!
 ) else (
        echo Fail.
	pause>nul
	exit
 )
plink -ssh -pw %password% -P %port% %user%@%sever% -batch "ls -l %remotePath%"

set /p "input=Test 22.3GB demo or not?(y/n)"
if /i %input% == "y" (
	start cmd /k plink -ssh -pw %password% -P %port% %user%@%sever%  -batch "top -1"
	plink -ssh -pw %password% -P %port% %user%@%sever%  -batch "%remotePath%/run_demo.sh"
	if %errorlevel% equ 0 (
        	echo Succeed!
   	 ) else (
        	echo Fail.
		pause>nul
		exit
    	)
	plink -ssh -pw %password% -P %port% %user%@%sever%  -batch "cat %remotePath%/ref_domo_result.txt"
)

endlocal

pause>nul