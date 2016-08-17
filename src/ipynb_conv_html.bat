@ECHO off
TITLE IPython Notebook [%DATE%] [%TIME%]
COLOR 8b

SET folder_bse=%~dp0
SET folder_src="%folder_bse%notebooks"
SET folder_dst="%folder_src%\..\nbconvert_html"

IF EXIST "%folder_dst%" RMDIR "%folder_dst%" /S/Q
MD "%folder_dst%"
CD /D "%folder_dst%"

jupyter-nbconvert --to=html "%folder_src%\*.ipynb"

PAUSE