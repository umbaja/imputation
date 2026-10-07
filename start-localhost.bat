@echo off
REM start-localhost.bat — dvojklik spusti Genome Converter appku cez WSL a otvori prehliadac.
setlocal
cd /d "%~dp0"

where wsl >nul 2>nul || (
  echo WSL sa nenasiel. Nainstaluj ho prikazom: wsl --install
  pause
  exit /b 1
)

echo ^>^> Startujem server vo WSL na http://127.0.0.1:8000
echo ^>^> Toto okno nezatvaraj — Ctrl+C ukonci server.
echo.

REM po 5 sekundach otvor prehliadac (server medzitym nabehne)
start "" /b cmd /c "timeout /t 5 >nul & start http://127.0.0.1:8000"

wsl -e bash -lc "cd \"$(wslpath -a '%~dp0')\" && bash serve.sh"

echo.
echo ^>^> Server sa ukoncil.
pause
