@echo off
setlocal enabledelayedexpansion

echo [MolGate] Stopping listeners on ports 5000, 5050, 8767...

for %%P in (5000 5050 8767) do (
  for /f "tokens=5" %%I in ('netstat -ano ^| findstr /R /C:":%%P .*LISTENING"') do (
    echo Killing PID %%I on port %%P
    taskkill /PID %%I /F >nul 2>nul
  )
)

echo [Done]
endlocal

