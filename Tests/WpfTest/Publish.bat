@echo off
setlocal
cd /d "%~dp0"

call Publish-small.bat
if errorlevel 1 goto error

call Publish-full-x64.bat
if errorlevel 1 goto error

call Publish-full-x32.bat
if errorlevel 1 goto error

pause
goto :eof

:error
chcp 65001 >nul
echo Публикация прервана из-за ошибки.
pause
exit /b 1
