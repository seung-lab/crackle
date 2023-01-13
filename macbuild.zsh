#!/bin/zsh

source ~/.zprofile

function compile_wheel {
	workon crackle$1
	pip install oldest-supported-numpy -r requirements.txt
	python setup.py bdist_wheel
}

compile_wheel 38
compile_wheel 39
compile_wheel 310
compile_wheel 311