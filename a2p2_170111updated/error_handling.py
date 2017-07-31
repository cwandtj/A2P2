__author__ = 'ones'
import os, signal



errors = {'Ones':['.FALSE.', '.TRUE.']}


def tail(f, n, offset=0):
  stdin,stdout = os.popen2("tail -n "+str(n)+ str(offset)+" "+f)
  stdin.close()
  lines = stdout.readlines()
  stdout.close()
  return lines


def check_error_logfile(directory, pid):
    incar_f = open(directory+'/INCAR', 'r')
    incar_lines = incar_f.readlines()

    log_f = directory+'/log'
    lines_new = tail(log_f, 5, 0)
    error_found = False
    lines_old = []

    while (lines_new != lines_old):

        for line in lines_new:
            for key in errors.keys():
                count = 0
                if (key in line):
                    error_found = True

                    for incar_line in incar_lines:
                        if 'LWAVE' in incar_line:
                            changed_line = ''
                            l = incar_line.split("=")
                            if l[1].strip() == errors[key][0]:
                                l[1] = errors[key][1]
                            elif l[1].strip() == errors[key][1]:
                                l[1] = errors[key][0]

                            changed_line += l[0] + ' = ' +l[1]

                        incar_lines[count] = changed_line
                        count += 1



        lines_old = lines_new
        lines_new = tail(log_f, 5, 0)

    changed_incar_f = open(directory+'/INCAR', 'w')
    changed_incar_f.writelines(incar_lines)

    if error_found:
        os.kill(pid, signal.SIGQUIT)
        os.chdir(directory)
        runjob = 'qsub runjob.sh'
        os.system(runjob)


