@echo off

for %%d in (bin obj .vs _ReSharper.Caches) do for /f %%f in ('dir /s /b /d /a %%d') do rd /s /q "%%f"

pause