@echo off
call "C:\Program Files\Microsoft Visual Studio 9.0\Common7\Tools\vsvars32.bat"

g77 -c *.f -fno-second-underscore -fno-leading-underscore
cl -c *.cpp
mkdir obj.win32
move /Y *.o obj.win32\.
move /Y *.obj obj.win32\.
del /F obj.win32\libXSModel.lib
ar -r obj.win32\libXSModel.lib obj.win32\*.o obj.win32\*.obj
