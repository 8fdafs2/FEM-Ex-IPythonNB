@echo off
TITLE IPython Notebook [%DATE%] [%TIME%]
COLOR 8b
CD /D %~dp0\notebooks
jupyter-nbconvert --execute --inplace --ExecutePreprocessor.timeout=600 *.ipynb
PAUSE