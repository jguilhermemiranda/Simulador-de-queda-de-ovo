@echo off
REM Gera egg_drop_lab_v3.exe na pasta do projeto (rodar no Windows).
python -m pip install pyinstaller numpy matplotlib pillow
if errorlevel 1 goto erro

python -m PyInstaller --noconfirm --clean --onefile --windowed ^
  --name egg_drop_lab_v3 ^
  --distpath . ^
  --workpath .build-work ^
  --specpath .build-spec ^
  --icon "%CD%\marca_idy_guy.ico" ^
  --add-data "%CD%\assets;assets" ^
  egg_drop_lab.py
if errorlevel 1 goto erro

echo.
echo OK: egg_drop_lab_v3.exe
exit /b 0
:erro
echo.
echo FALHOU. Veja a mensagem acima.
exit /b 1
