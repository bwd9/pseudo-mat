from random import choice, random
import os

current_version = 11
run_path = os.environ["HOME"]+'/mat'+str(current_version)
#match run_path name to main...

if not os.path.exists(run_path):
	#mkmatdir = 'mkdir '+ matpath
	os.popen('mkdir ' + run_path, 'w', 1)

test = open(os.path.join(os.path.abspath(mat_path), 'testfile', 'a')
test.write('write test')
