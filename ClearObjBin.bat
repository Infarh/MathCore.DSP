@echo off

for /f %%f in ('dir /s /b obj') do (del /f /s /q %%f > nul)
for /f %%f in ('dir /s /b bin') do (del /f /s /q %%f > nul)

for /f %%f in ('dir /s /b obj') do (rmdir /s /q %%f > nul)
for /f %%f in ('dir /s /b bin') do (rmdir /s /q %%f > nul)