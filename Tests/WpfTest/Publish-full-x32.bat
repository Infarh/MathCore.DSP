@echo off
setlocal
cd /d "%~dp0"

dotnet publish -v d /p:PublishProfile=Properties/PublishProfiles/Full-x32.pubxml

if errorlevel 1 goto error

set PublishDir=!release\WpfTest-Full-x32
set ArchiveName=!release\WpfTest-Full-x32.zip

if exist "%ArchiveName%" del "%ArchiveName%"

powershell -NoProfile -Command "Compress-Archive -Path '%PublishDir%\*' -DestinationPath '%ArchiveName%' -Force"

chcp 65001 >nul

echo Публикация завершена. Архив создан: %ArchiveName%
goto :eof

:error
set ExitCode=%errorlevel%
chcp 65001 >nul
echo Ошибка публикации. Код завершения: %ExitCode%
exit /b %ExitCode%
