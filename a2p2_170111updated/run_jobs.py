__author__ = 'ones'

import os, sys, threading
import error_handling

working_dir = os.getcwd()
directories = []
zombies = []

if len(sys.argv) == 1:
    os.system('find -mindepth 3 -type d > compounds_directories')


elif len(sys.argv) == 2:
    root = sys.argv[1]
    depth = root.count('/')
    os.system('find ./'+root+' -mindepth '+str(2-depth)+' -type d > compounds_directories')
else:
    exit("No more than 1 argument")


dir_f = open('compounds_directories', 'r')
lines = dir_f.readlines()
dir_f.close()
max_depth = 0
for line in lines:
    depth = len(filter(None, line.split('/')))
    if depth > max_depth:
        max_depth = depth

for line in lines:
    line = line.strip()[1:]
    if (len(filter(None, line.split('/'))) == max_depth-1):
        directories.append(working_dir+line)


for directory in directories:
    os.chdir(directory)
    runjob = 'qsub runjob.sh'
    os.system(runjob)
    

    # t1 = threading.Thread(target= runjob, args=(directory, ))
    # t1.setDaemon(True)
    # t1.start()

    # t2 = threading.Thread(target= handle_errors, args=(directory, ))
    # t2.setDaemon(True)
    # t2.start()





#
#     process = os.fork()
#     if process == 0:
#         os.chdir(directory)
#         runjob = 'qsub runjob.sh'
#         os.system(runjob)
#
#
#     else:
#         error_handling.check_error_logfile(directory, process)
#
#
#     zombies.append(process)
#
# for zombie in zombies:
#     os.waitpid(zombie, 0)












